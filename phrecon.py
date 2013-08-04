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

import argparse
import os
import sys
import re
import sqlite3
import Bio
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import Bio.Alphabet

VERSION = "3.0"

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
        logfile.write(msg)

def print_all(msg):
    '''Prints output to STDOUT and logfile in one call.
       If the quiet funnction is set, then only log to file.

       Parameters:
       msg = Text to print/log.
    '''
    if not quiet:
        print msg
    logfile.write(msg)

def write_fasta(strainid,sequenceString):
    '''Writes out the FASTA formatted entry for the strainid
       with the given sequence String to the outfile
    '''

    outfile.write("\n")
    outfile.write(">{}\n".format(strainid))
    # Chop up the sequence to 60 characters per line.
    fastaRow = [sequenceString[i:i+60] for i in range(0,len(sequenceString),60) ]
    for fastaLine in fastaRow:
        outfile.write("{}\n".format(fastaLine))
    




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
parser.add_argument('-dbfile','--dbfile', type=str, nargs=1, help="Working database name (default is to load into memory)")
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

print_debug("Detected Biopython version {}\n".format(Bio.__version__))

# Setup in-memory (or disk, if --dbfile is used) database.
if args.dbfile != None:
    dbfilename = args.dbfile[0]
    if os.path.isfile(dbfilename):
        print "*** ERROR database file {0} already exists.  Please remove or use a different filename.\n".format(dbfilename)
        sys.exit(1)
    else:
        dbconn = sqlite3.connect(dbfilename)
        dbconn.execute("PRAGMA synchronous = OFF")
        dbconn.execute("PRAGMA temp_store = MEMORY")
        dbconn.execute("PRAGMA journal_mode = MEMORY")
else:
    dbconn = sqlite3.connect(':memory:')

dbcursor = dbconn.cursor()

# Create the table for the ref file data (COLUMNS: locus, base)
dbcursor.execute('''CREATE TABLE REF_DATA (locus integer PRIMARY KEY ASC, base text NOT NULL)''')
dbcursor.execute('''CREATE INDEX REF_DATA_LOCUS_BASE_IDX ON REF_DATA (locus,base)''')

# Create the table for the snp file data (COLUMNS: strainid, locus, base)
dbcursor.execute('''CREATE TABLE SNP_DATA (strainid text, locus integer, base text NOT NULL, PRIMARY KEY(strainid,locus))''')

# Create indexes.
dbcursor.execute('''CREATE INDEX SNP_DATA_LOCUS_IDX ON SNP_DATA (locus)''')
dbcursor.execute('''CREATE INDEX SNP_DATA_STRAINID_IDX ON SNP_DATA (strainid)''')
dbcursor.execute('''CREATE INDEX SNP_DATA_STRAINID_LOCUS_BASE_IDX ON SNP_DATA (strainid,locus,base)''')

dbconn.commit()


print_all("Loading loci entries from reference file...\n")
# Read in and store reference data in REF_DATA table.
lineNumber=0

# Expect reference data to be in a two-column format: Loci Base
refRe = re.compile("^(?P<locus>\d+)\s+(?P<base>[A-Z]+)$")
for line in reffile:
    lineNumber += 1


    refMatch = refRe.match(line)
    if refMatch:
        locus = refMatch.group('locus')
        base = refMatch.group('base')
        # Check for conflict.  Fail if mismatched base at locus, ignore if identical base at locus, insert if no row exists.
        dbcursor.execute("SELECT base FROM REF_DATA where locus = ?",(locus,))
        dbresult = dbcursor.fetchone()
        if dbresult and dbresult[0] != base:
            print_all("*** ERROR: Duplicate loci found! Read base {} at loci {}, which collides with previous base {} at same loci!  On line {} of reference file.\n".format(base,locus,dbresult[0],lineNumber))
            print "Failed."
            sys.exit(1)
        elif dbresult and dbresult[0] == base:
            print_debug("Skipping duplicate (but identical) locus/base at line {} for locus {} and base {}\n".format(lineNumber,locus,base))
        else:
            dbconn.execute("INSERT INTO REF_DATA (locus, base) VALUES (?,?)",(locus,base))
            dbconn.commit()
    else:
        print_all("*** ERROR: Reference file format not recognized.  Expected two columns -- [Loci (numeric)] [Base (A-Z)].  Got {} at line {} instead.\n".format(line,lineNumber))
        print "Failed.\n"
        sys.exit(1)

print
reffile.close()
print_all("Read {} loci entries from reference file.\n".format(lineNumber))

print_all("Loading loci entries from strain snp data file...\n")
# Read in and store SNP data in SNP_DATA table.
lineNumber=0
warnings=0
strainCount = 0

# Expect reference data to be in a three-column format: StrainID Loci Base
snpRe = re.compile("^(?P<strainid>\S+)\t(?P<locus>\d+)\t(?P<base>[A-Z]+)$")
with dbconn:
    currentStrain = ""
    for line in infile:
        lineNumber += 1

        snpMatch = snpRe.match(line)

        # Quit on bad lines.  Bad data!
        if snpMatch == None:
            print_all("*** ERROR: Bad data format on line {} of SNP loci file.  Expect three columns [Strain ID (Alphanumeric)] [Loci (Numeric)] [Base (alphabetic)]\nGot {} instead.\n".format(lineNumber,line))
            print "Failed.\n"
            sys.exit(1)

        strainid = snpMatch.group('strainid')
        locus = snpMatch.group('locus')
        base = snpMatch.group('base')

        if currentStrain != strainid:
            strainCount += 1
            currentStrain = strainid

        # Sanity checks - These just generate warnings
        # Check that the locus exists in the reference base.
        dbcursor.execute("SELECT * from REF_DATA where locus = ?",(locus,))
        dbresult = dbcursor.fetchone()
        if not dbresult:
            print_all("*** WARNING: Encountered SNP loci {}, which does NOT exist in base ref sequence! (line {} of SNP input file) Skipping...\n".format(locus,lineNumber))
            warnings += 1
            continue

        # Check that the locus doesn't already exist for this strain
        dbcursor.execute("SELECT * from SNP_DATA where strainid = ? and locus = ?",(strainid,locus,))
        dbresult = dbcursor.fetchone()
        if dbresult:
            print_all("*** WARNING: Duplicate loci found for same strain id! (line {} of SNP input file, colliding on loci {}) Skipping (will use existing loci)...\n".format(lineNumber,locus))
            warnings += 1
            continue


        # Load into database
        print_debug("Diag: read in {}, {}, {}\n".format(strainid,locus,base))
        dbconn.execute("INSERT INTO SNP_DATA (strainid, locus, base) VALUES (?,?,?)",(strainid,locus,base))

print
infile.close()
print_all("Read {} loci entries from snp data file.\n".format(lineNumber))

seqRecordList = []

# Store the reference base sequence as the first output in the sequence record array
print_all("Generating reference sequence record.\n")
dbcursor.execute('''select group_concat(base,"") from REF_DATA order by locus asc;''')
dbresult = dbcursor.fetchone()
if dbresult:
    refSequenceString = dbresult[0]
    refSeq = Bio.Seq.Seq(refSequenceString,Bio.Alphabet.Alphabet)
    refSeqRecord = Bio.SeqRecord.SeqRecord(refSeq,id="REF",description="")
    seqRecordList.append(refSeqRecord)
else:
    print_all("*** ERROR: Failed to generate reference sequence record.  No data?\n")
    sys.exit(1)
    
# Reconstruct the input strains and store those into the sequence record array
straincursor = dbconn.cursor()
straincursor.execute('''select distinct strainid from SNP_DATA order by strainid;''')
for strainRow in straincursor.fetchall():
    currentStrainID = strainRow[0]
    print_all("Now reconstructing strain {}".format(currentStrainID))

    dbcursor.execute('''select group_concat(coalesce(b.base,a.base),"") from REF_DATA a LEFT OUTER JOIN SNP_DATA b on a.locus = b.locus and b.strainid=?  order by a.locus asc;''',(currentStrainID,))
    dbresult = dbcursor.fetchone()

    snpSequenceString = dbresult[0]
    snpSeq = Bio.Seq.Seq(snpSequenceString,Bio.Alphabet.Alphabet)
    snpSeqRecord = Bio.SeqRecord.SeqRecord(snpSeq,id=currentStrainID,description="")
    seqRecordList.append(snpSeqRecord)

# Write it all out
print_all("Writing output file...\n")
Bio.SeqIO.write(seqRecordList,outfilename,"fasta")

print_all("\n{} strains processed from SNP loci input file.\n".format(strainCount));
print_all("\n{} warnings were detected.\n".format(warnings));
print "Log file for this run can be found in {}\n".format(logfilename)


