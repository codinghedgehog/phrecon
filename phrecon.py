#!/usr/bin/python
#
# Phrecon (Phylo Reconstructor)
#
# This script reconstructs the base sequence using an SNP loci and substituting its bases in at the given locations in the
# provided reference base sequence.
#
# The base reference should only have TWO columns (and no headers), consisting of just the space-delimited locus and base columns.
#
# The SNP input file may contain multiple SNP "fragments" (ids), and phrecon will generate a new sequence for each one
# and write them out in FASTA format into the same output file.
#
# SNP file input should contain three columns and no headers, space-delimited:
# SNP_ID (i.e. query id) | Base position | Base
#
# So something like:
# A12 1045  G
# A12 4056  A
# A12 13004 T
# A35 4 A
# A35 401 C
#
# The above example contains two query SNPs, so phrecon will generate a FASTA file containing the full sequence for each.
#
# Prephix's SNP and base ref output files can be used by phrecon as input.
#
# If the SNP line has a locus of -1 and base "-", and it is the ONLY line for that strainid, then it is treated as
# an empty SNP -- i.e. the same sequence as the ref base (since there are no SNPs reported for that strain).
#
# Usage: phrecon.py <reference base file> <SNP loci input file> [-debug] [-quiet]
#
# 8/2/2013 - Version 3 - Andrew Pann, Ported version 2.1 of Perl-based Phrecon to Python, using Biopython and Sqlite3.
# 8/4/2013 - Version 3.1 - Andrew Pann, Reverted to normal port using dictionaries and Biopython (no database).
# 8/4/2013 - Version 4.0 - Andrew Pann, Added use of multiprocessing.
# 8/2/2013 - Version 4.1 - Andrew Pann, Added handling of empty strain SNP data.
# 4/13/2014 - Version 4.2 - Andrew Pann, Fixed handling of single locus SNP data inputs.

import argparse
import os
import sys
import re
import Bio
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet
import multiprocessing
import pprint
import signal

VERSION = "4.2.0"

#####################
# UTILITY FUNCTIONS
#####################

def print_debug(msg):
    '''Prints output to STDOUT and logfile if debug flag is set.
       If the quiet funnction is also set, then only log to file.

       Parameters:
       msg = Text to print/log.
    '''
    if debug:
        if not quiet:
            print msg
        logfile.write("{}\n".format(msg))

def print_all(msg):
    '''Prints output to STDOUT and logfile in one call.
       If the quiet funnction is set, then only log to file.

       Parameters:
       msg = Text to print/log.
    '''
    if not quiet:
        print msg
    logfile.write("{}\n".format(msg))


########################################
# FUNCTIONS RELATED TO MULTIPROCESSING
########################################

def mp_reconstruct_strain(myStrainID,myStrainData,myRefData):
    '''This function is used by the multiprocessing feature.
       It takes in three arguments: strainid name, strain id's
       dictionary object (locus=>base), the reference base 
       dictionary data (locus=>base).

       This function returns a two-tuple of the strainid name and
       the string representing the reconstructed snp sequence,
       or an Exception object as the second two-tuple value if
       there is an error.

       If the strain data has one entry of locus -1 and base '-', 
       then it is assumed to be from an empty SNP file (no SNPs),
       so will return the base reference sequence.
    '''
    try:
        if len(myStrainData) == 1 and 0 in myStrainData.keys() and myStrainData[0] == "-1":
            strainData = myRefData.copy()
        else:
            strainData = myRefData.copy()
            strainData.update(myStrainData)

        strainKeys = strainData.keys()
        strainKeys.sort()
        snpSequenceString = "".join([strainData[locus] for locus in strainKeys])
    except Exception as error:
        return (myStrainID,error)
    else:
        return (myStrainID,snpSequenceString)

def mp_init_worker():
    '''Initializer function for Pool. Called by each worker process, and will
       make them ignore interrupts, to be handled by the main thread.
    '''
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def mp_store_result_callback(result):
    '''Callback function for Pool's async calls.  Stores the result
       in the result table.
    '''
    logFileLock.acquire()
    if isinstance(result[1],Exception):
        print_all("*** ERROR: Exception thrown while processing strain {}: {}".format(result[0],pprint.pformat(result[1])))
    else:
        print_all("Now reconstructing strain {}".format(result[0]))
        resultTable[result[0]] = result[1]
    logFileLock.release()
        

########################
# MAIN
########################


if __name__ == '__main__':

    print "\nPhrecon (Phylo Reconstructor) v{}\n\n".format(VERSION)

    # Setup argument parsing

    parser = argparse.ArgumentParser()
    parser.add_argument('reffilename', help="Reference base file name with line format [LOCUS] [BASE]")
    parser.add_argument('snpfilename', help="SNP data file name with line format [STRAIN_ID] [LOCUS] [BASE]")
    parser.add_argument('-debug','--debug', action="store_true", help="Enable debugging output")
    parser.add_argument('-quiet','--quiet', action="store_true", help="Run in quiet mode (suppress most screen output)")
    parser.add_argument('-cpus','--cpus', type=int, help="Indicates how many CPUs to use during multiprocessing.  Default is the detected number of CPUs from Python's multiprocessing module.  Setting this to 1 basically disables multiprocessing.",default=multiprocessing.cpu_count())
    args = parser.parse_args()

    reffile = open(args.reffilename,"r")
    infile = open(args.snpfilename,"r")
    logfilename = "{}.reconstructed.log".format(args.snpfilename)
    logfile = open(logfilename,"w")

    #outfile is not needed since BioPython will write to it.  Just need to generate the filename here.
    outfilename = "{}.reconstructed".format(args.snpfilename)

    debug = args.debug
    quiet = args.quiet

    cpu_num=args.cpus

    if debug:
        print "Producing debug output.\n"

    if quiet:
        print "Producing quiet (no stdout) output.  Logfile is still generated.\n"

    print_all("CPUs to use: {}".format(cpu_num))

    refData = {} # This is a dictionary of the ref base. Key is locus, value is base.

    print_all("Loading loci entries from reference file...")
    # Read in and store reference data in REF_DATA table.
    lineNumber=0

    # Expect reference data to be in a two-column format: Loci Base
    refRe = re.compile("^(?P<locus>\d+)\s+(?P<base>[A-Z]+)$")
    for line in reffile:
        lineNumber += 1


        refMatch = refRe.match(line)
        if refMatch:
            locus = int(refMatch.group('locus'))
            base = refMatch.group('base')
            # Check for conflict.  Fail if mismatched base at locus, ignore if identical base at locus, insert if no row exists.
            if locus in refData:
                if refData[locus] != base:
                    print_all("*** ERROR: Duplicate loci found! Read base {} at loci {}, which collides with previous base {} at same loci!  On line {} of reference file.".format(base,locus,refData[locus],lineNumber))
                    print "Failed."
                    sys.exit(1)
                else:
                    print_debug("Skipping duplicate (but identical) locus/base at line {} for locus {} and base {}\n".format(lineNumber,locus,base))
            else:
                refData[locus] = base
        else:
            print_all("*** ERROR: Reference file format not recognized.  Expected two columns -- [Loci (numeric)] [Base (A-Z)].  Got {} at line {} instead.".format(line,lineNumber))
            print "Failed.\n"
            sys.exit(1)

    print
    reffile.close()
    print_all("Read {} loci entries from reference file.".format(lineNumber))

    print_all("Loading loci entries from strain snp data file...")
    # Read in SNP data and store in the snpData dictionary
    # The snpData dictionary has strainIDs as the key, and the values are
    # a dictionary of locus keys and base values for that strain.
    lineNumber=0
    warnings=0
    strainCount = 0
    strainNameList = [] # Store strain IDs in a list to preserve order when writing to file later.
    # Expect snp data file lines to be in a three-column format: StrainID Loci Base
    snpRe = re.compile("^(?P<strainid>\S+)\t(?P<locus>\d+|(-1))\t(?P<base>[A-Z]+|-)$")
    currentStrain = ""
    snpData = {}
    for line in infile:
        lineNumber += 1

        snpMatch = snpRe.match(line)

        # Quit on bad lines.  Bad data!
        if snpMatch == None:
            print_all("*** ERROR: Bad data format on line {} of SNP loci file.  Expect three columns [Strain ID (Alphanumeric)] [Loci (Numeric)] [Base (alphabetic)]\nGot {} instead.".format(lineNumber,line))
            print "Failed.\n"
            sys.exit(1)

        strainid = snpMatch.group('strainid')
        locus = int(snpMatch.group('locus'))
        base = snpMatch.group('base')

        if currentStrain != strainid:
            strainCount += 1
            currentStrain = strainid
            strainNameList.append(currentStrain)
            snpData[currentStrain] = {}

        # Sanity checks - These just generate warnings
        # Check that the locus exists in the reference base.
        if not locus in refData:
            print_all("*** WARNING: Encountered SNP loci {}, which does NOT exist in base ref sequence! (line {} of SNP input file) Skipping...".format(locus,lineNumber))
            warnings += 1
            continue

        # Check that the locus doesn't already exist for this strain
        if locus in snpData[currentStrain]:
            print_all("*** WARNING: Duplicate loci found for same strain id! (line {} of SNP input file, colliding on loci {}) Skipping (will use existing loci)...".format(lineNumber,locus))
            warnings += 1
            continue

        # Load into SNP dictionary
        print_debug("Diag: read in {}, {}, {}\n".format(strainid,locus,base))
        snpData[currentStrain][locus]=base

    print
    infile.close()
    print_all("Read {} loci entries from snp data file.".format(lineNumber))

    seqRecordList = []

    # Store the reference base sequence as the first output in the sequence record array
    print_all("Generating reference sequence record.")

    refKeys = refData.keys()
    refKeys.sort()
    refSequenceString = "".join([ refData[locus] for locus in refKeys])
    refSeq = Bio.Seq.Seq(refSequenceString,Bio.Alphabet.Alphabet)
    refSeqRecord = Bio.SeqRecord.SeqRecord(refSeq,id="REF",description="")
    seqRecordList.append(refSeqRecord)
        
    ##################
    # MULTIPROCESSING -- Reconstruct the input strains in parallel and store those into the sequence record array.
    ##################

    # Setup a lock -- this is just to neatly interleave output to the log file (otherwise missing or garbles lines can occur).
    logFileLock = multiprocessing.Lock()
    resultTable = {} # Dictionary of strainID => base sequence string.  Filled by multiprocesing call.
    mpPool = multiprocessing.Pool(processes=cpu_num,initializer=mp_init_worker)
    print_all("Submitting strain reconstruction requests to pool...")
    try:
        for strainID in strainNameList: 
            mpPool.apply_async(mp_reconstruct_strain,args=(strainID,snpData[strainID],refData),callback=mp_store_result_callback)

    except KeyboardInterrupt:
        print "Aborting..."
        mpPool.terminate()
        mpPool.join()
    else:
        print_all("Waiting for working processes to return...")
        mpPool.close()
        mpPool.join()

    # Process results from multiprocessing.
    # Namely this means generating the list of SeqRecords (in same order as the snp input file)
    # to write to the output file with Biopython.
    print_all("Processing results...")
    for strainID in strainNameList:
        snpSequenceString = resultTable[strainID]
        snpSeq = Bio.Seq.Seq(snpSequenceString,Bio.Alphabet.Alphabet)
        snpSeqRecord = Bio.SeqRecord.SeqRecord(snpSeq,id=strainID,description="")
        seqRecordList.append(snpSeqRecord)

    # Write it all out.
    print_all("Writing output file...")
    Bio.SeqIO.write(seqRecordList,outfilename,"fasta")

    print_all("\n{} strains processed from SNP loci input file.".format(strainCount));
    print_all("\n{} warnings were detected.".format(warnings));
    print "Log file for this run can be found in {}\n".format(logfilename)


