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
# Usage: phrecon.py <reference base file> <SNP loci input file> [-debug] [-quiet]
#
# 8/2/2013 - Version 3 - Andrew Pann, Ported version 2.1 of Perl-based Phrecon to Python, using Biopython and Sqlite3.
# 8/4/2013 - Version 3.1 - Andrew Pann, Added multiprocessing to leverage multi-CPU performance.

import argparse
import os
import sys
import re
import Bio
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet

VERSION = "3.1"

#############
# FUNCTIONS
#############

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


########################
# MAIN
########################

print "\nPhrecon (Phylo Reconstructor) v{}\n\n".format(VERSION)

# Setup argument parsing

parser = argparse.ArgumentParser()
parser.add_argument('reffilename', help="Reference base file name with line format [LOCUS] [BASE]")
parser.add_argument('snpfilename', help="SNP data file name with line format [STRAIN_ID] [LOCUS] [BASE]")
parser.add_argument('-debug','--debug', action="store_true", help="Enable debugging output")
parser.add_argument('-quiet','--quiet', action="store_true", help="Run in quiet mode (suppress most screen output)")
args = parser.parse_args()

reffile = open(args.reffilename,"r")
infile = open(args.snpfilename,"r")
logfilename = "{}.reconstructed.log".format(args.snpfilename)
logfile = open(logfilename,"w")

#outfile is not needed since BioPython will write to it.  Just need to generate the filename here.
outfilename = "{}.reconstructed".format(args.snpfilename)

debug = args.debug
quiet = args.quiet

if debug:
    print "Producing debug output.\n"

if quiet:
    print "Producing quiet (no stdout) output.  Logfile is still generated.\n"


refData = {} # This is a dictionary of the ref base. Key is locus, value is base.

print_all("Loading loci entries from reference file...\n")
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
                print_all("*** ERROR: Duplicate loci found! Read base {} at loci {}, which collides with previous base {} at same loci!  On line {} of reference file.\n".format(base,locus,refData[locus],lineNumber))
                print "Failed."
                sys.exit(1)
            else:
                print_debug("Skipping duplicate (but identical) locus/base at line {} for locus {} and base {}\n".format(lineNumber,locus,base))
        else:
            refData[locus] = base
    else:
        print_all("*** ERROR: Reference file format not recognized.  Expected two columns -- [Loci (numeric)] [Base (A-Z)].  Got {} at line {} instead.\n".format(line,lineNumber))
        print "Failed.\n"
        sys.exit(1)

print
reffile.close()
print_all("Read {} loci entries from reference file.\n".format(lineNumber))

print_all("Loading loci entries from strain snp data file...\n")
# Read in SNP data and store in the snpData dictionary
# The snpData dictionary has strainIDs as the key, and the values are
# a dictionary of locus keys and base values for that strain.
lineNumber=0
warnings=0
strainCount = 0
strainNameList = [] # Store strain IDs in a list to preserve order when writing to file later.
# Expect reference data to be in a three-column format: StrainID Loci Base
snpRe = re.compile("^(?P<strainid>\S+)\t(?P<locus>\d+)\t(?P<base>[A-Z]+)$")
currentStrain = ""
snpData = {}
for line in infile:
    lineNumber += 1

    snpMatch = snpRe.match(line)

    # Quit on bad lines.  Bad data!
    if snpMatch == None:
        print_all("*** ERROR: Bad data format on line {} of SNP loci file.  Expect three columns [Strain ID (Alphanumeric)] [Loci (Numeric)] [Base (alphabetic)]\nGot {} instead.\n".format(lineNumber,line))
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
        print_all("*** WARNING: Encountered SNP loci {}, which does NOT exist in base ref sequence! (line {} of SNP input file) Skipping...\n".format(locus,lineNumber))
        warnings += 1
        continue

    # Check that the locus doesn't already exist for this strain
    if locus in snpData[currentStrain]:
        print_all("*** WARNING: Duplicate loci found for same strain id! (line {} of SNP input file, colliding on loci {}) Skipping (will use existing loci)...\n".format(lineNumber,locus))
        warnings += 1
        continue

    # Load into SNP dictionary
    print_debug("Diag: read in {}, {}, {}\n".format(strainid,locus,base))
    snpData[currentStrain][locus]=base

print
infile.close()
print_all("Read {} loci entries from snp data file.\n".format(lineNumber))

seqRecordList = []

# Store the reference base sequence as the first output in the sequence record array
print_all("Generating reference sequence record.\n")

refKeys = refData.keys()
refKeys.sort()
refSequenceString = "".join([ refData[locus] for locus in refKeys])
refSeq = Bio.Seq.Seq(refSequenceString,Bio.Alphabet.Alphabet)
refSeqRecord = Bio.SeqRecord.SeqRecord(refSeq,id="REF",description="")
seqRecordList.append(refSeqRecord)
    
# Reconstruct the input strains and store those into the sequence record array
for strainID in strainNameList: 
    print_all("Now reconstructing strain {}".format(strainID))
    strainData = refData.copy()
    strainData.update(snpData[strainID])

    strainKeys = strainData.keys()
    strainKeys.sort()
    snpSequenceString = "".join([strainData[locus] for locus in strainKeys])
    snpSeq = Bio.Seq.Seq(snpSequenceString,Bio.Alphabet.Alphabet)
    snpSeqRecord = Bio.SeqRecord.SeqRecord(snpSeq,id=strainID,description="")
    seqRecordList.append(snpSeqRecord)

# Write it all out
print_all("Writing output file...\n")
Bio.SeqIO.write(seqRecordList,outfilename,"fasta")

print_all("\n{} strains processed from SNP loci input file.\n".format(strainCount));
print_all("\n{} warnings were detected.\n".format(warnings));
print "Log file for this run can be found in {}\n".format(logfilename)


