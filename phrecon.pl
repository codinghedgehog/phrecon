#!/usr/bin/perl
#
# Phrecon (Phylo Reconstructor)
#
# This script reconstructs the base sequence using an SNP loci and substituting its bases in at the given locations in the 
# provided reference base sequence.
#
# The base reference should only have TWO columns (and no headers), consisting of just the base location and base.
#
# The SNP input file may contain multiple SNP "fragments" (ids), and phrecon will generate a new sequence for each one
# and write them out in FASTA format into the same output file.
#
# SNP file input should contain three columns and no headers:
# SNP_ID (i.e. query id) | Base position | Base
#
# So something like:
# A12	1045	G
# A12	4056	A
# A12 13004	T
# A35	4	A
# A35	401	C
#
# The above example contains two query SNPs, so phrecon will generate a FASTA file containing the full sequence for each.
#
# Prephix's SNP and base ref output files can be used by phrecon as input.
#
# Usage: phrecon.pl <reference base file> <SNP loci input file>
#
# 4/26/2011 - Andrew Pann & S. Wesley Long- Initial development.
# 1/26/2012 - Andrew Pann, made errors of duplicate loci and nonexistent ref loci to be warnings, instead of hard errors.  Also allow alphanumeric strain IDs.
# 1/27/2012 - Andrew Pann, added log file and support for -debug and -quiet flag.
# 2/02/2012 - Andrew Pann, added failure message on errors, even in quiet mode.  No silent failures!

use strict;

my $VERSION="2.1";

print "\nPhrecon (Phylo Reconstructor) v$VERSION\n\n";

if ($#ARGV < 1){
  print "Usage: $0 <reference base file> <SNP loci input file> [-debug] [-quiet]\n";
  exit 1;
}

my $infile;
my $reffile;
my $outfile;
my $logfile;
my %baseRefTable;
my %outputTable;
my %snpTable;
my $i=0;
my $currentStrain=0;
my $loci;
my $strainCount=0;
my $charCount=0;
my $warnings=0;
my $debug="N";
my $quiet="N";
my $arg_num=0;

open($reffile,"<",$ARGV[0]) or die "Unable to open reference file $ARGV[0] for reading!  $!\n";
open($infile,"<",$ARGV[1]) or die "Unable to open input file $ARGV[1] for reading!  $!\n";
open($outfile,">","$ARGV[1].reconstructed") or die "unable to open output file $ARGV[1].reconstructed for writing! $!\n";
open($logfile,">","$ARGV[1].reconstructed.log") or die "unable to open log file $ARGV[1].reconstructed.log for writing! $!\n";

# Process optional parameters
if ($#ARGV > 1){
  $arg_num=2;
  while ($arg_num <= $#ARGV){
    # Using if-else logic because given/when isn't compatible with< 5.10 and Switch module was removed from >5.14.
    if ($ARGV[$arg_num] eq "-debug"){
				$debug="Y";
        print "Producing debug output.\n";
    }
    elsif ($ARGV[$arg_num] eq "-quiet") { 
        print "Producing quiet (no stdout) output.  Log file is still generated.\n";
				$quiet="Y";
    }
    else{
        print "*** ERROR: Unknown command line argument $ARGV[$arg_num].\n";
        exit 1;
    }
    $arg_num++;
  }
}  

print_all("Loading loci entries from reference file...\n");
# Read in and store reference data in a hash table keyed off of loci number.
while (<$reffile>){
  $i++;
  #print "\rEntries read: $i";

  chomp;

  # Expect reference data to be in a two-column format: Loci Base
  if (/^([0-9]+)\s+([A-Z]+)$/){
    if (exists($baseRefTable{$1})){
      print_all("*** ERROR: Duplicate loci found! Read base $2 at loci $1, which collides with previous base $baseRefTable{$1} at same loci!  On line $i of reference file.\n");
			print "Failed.\n";
      exit 1;
    }
    else{
      $baseRefTable{$1}="$2";
    }
  }
  else{
    print_all("*** ERROR: Reference file format not recognized.  Expected two columns -- [Loci (numeric)] [Base (A-Z)].  Got \"$_\" at line $i instead.\n");
		print "Failed.\n";
    exit 1;
  }
}
print "\n";
close($reffile);
print_all("Read $i loci entries from reference file.\n");

# Force the reference table to be output as the first entry.
print_all("Writing reference sequence.\n");
$charCount = 0;
print $outfile ">REF\n";
foreach $loci (sort {$a <=> $b} keys %baseRefTable){

  print $outfile "$baseRefTable{$loci}";

  $charCount++;

  # Limit to 70 chars per line, for FASTA compliance.
  if ($charCount eq 70){
    print $outfile "\n";
    $charCount = 0;
  }

}
print $outfile "\n";

$i = 0;
# For each strain in the input file, output the full sequence.
# Method: Make a copy of the base reference hash table, and change its values as needed for each strain
# via the loci info from the SNP file.
while (<$infile>){
  $i++;

  chomp;

  if (!(/^(\S+)\t(\d+)\t([A-Z]+)$/)){
    print_all("*** ERROR: Bad data format on line $i of SNP loci file.  Expect three columns [Strain ID (Alphanumeric)] [Loci (Numeric)] [Base (alphabetic)]\nGot \"$_\" instead.\n");
		print "Failed.\n";
    exit 1;
  }

  my ($in_strain,$in_loci,$in_base) = split /\t/;

  print_debug("Diag: Read in $in_strain, $in_loci, $in_base\n");

  # Check if new strain encountered.   If so, dump out the sequence of the strain we just finished.
  # Then output new strain FASTA header and reset output hash table to start reconstructing the new strain.
  if ($in_strain ne $currentStrain){
    print_debug("Diag: New strain encountered: $in_strain\n");
    $strainCount++;
    $charCount = 0;

    # Dump out previous strain sequence that's been reconstructed up to this point.
    print_debug("Diag: Writing strain $currentStrain to output file.\n");
    foreach $loci (sort {$a <=> $b} keys %outputTable){

      print $outfile "$outputTable{$loci}";

      $charCount++;

      # Limit to 70 chars per line, for FASTA compliance.
      if ($charCount eq 70){
        print $outfile "\n";
	$charCount = 0;
      }

    }
    print $outfile "\n";

    # Set up for reconstruction of new strain.
    $currentStrain = $in_strain;
    print $outfile ">$currentStrain\n";
    %outputTable = %baseRefTable;
    print_all("Now reconstructing strain $currentStrain\n");

    %snpTable = ( );
  }

  if (!(exists($outputTable{$in_loci}))){
    print_all("*** WARNING: Encountered SNP loci $in_loci, which does NOT exist in base ref sequence! (line $i of SNP input file) Skipping...\n");
    $warnings++;
  }
  else{
    if (exists($snpTable{$in_loci})){
      print_all("*** WARNING: Duplicate loci found for same strain id! (line $i of SNP input file, colliding on loci $in_loci) Skipping (will use existing loci)...\n");
      $warnings++;
    }
    else{
      $outputTable{$in_loci} = $in_base;
      $snpTable{$in_loci} = $in_base;
    }
  }
}

# Dump out the last strain's sequence that's been reconstructed up to this point.
$charCount = 0;
foreach $loci (sort {$a <=> $b} keys %outputTable){

  print $outfile "$outputTable{$loci}";

  $charCount++;

  # Limit to 70 chars per line, for FASTA compliance.
  if ($charCount eq 70){
    print $outfile "\n";
    $charCount = 0;
  }

}
print $outfile "\n";

print_all("\n$strainCount strains processed from SNP loci input file.\n");
print_all("\n$warnings warnings were detected.\n");
print "Log file for this run can be found in $ARGV[1].reconstructed.log \n";

close($infile);
close($outfile);

#############
# FUNCTIONS
#############
sub print_debug(){
  # Prints output to STDOUT and log file if debug flag is set.  Otherwise nothing.
  # If the quiet function is also set, then only log to file.
  #
  # Parameters:
  # $1 = Text to print.
  #
  if ($debug eq "Y"){
    if ($quiet eq "N"){
      print "$_[0]";
    }
    print $logfile "$_[0]";
  }
}

sub print_all(){
  # Silly MUX'ed function to write output to both standard out and logfile in one call.
  # If the quiet function is set, then only log to file.
  #
  # Parameters:
  # $1 = Text to print.
  #
  if ($quiet eq "N"){
    print "$_[0]";
  }
  print $logfile "$_[0]";
}

