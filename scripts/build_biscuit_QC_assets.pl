#!/usr/bin/perl -w
use strict;
=head1
Author: Ian Beddows (ian.beddows at vai.org)
Date: 8/4/2020

Update Notes:

    25 Aug 2020 -
        - Updating header and usage
        - Count soft-masked CpGs (CG/Cg/cG/cg are all counted as CpGs)
        - Correct GC-content calculation
        - Clean up code for ease of maintaining
    21 Oct 2020 -
        - Update sorting for GC-content windows

Description:

This program takes an input FASTA file (reference genome) and creates 3 output
files in the provided output directory. The three files and a rough description
of how they are created are:

   1. cpg.bed.gz
        - Find positions of all CpGs in genome reference FASTA.
        - Put positions in BED format. You only need the chromosome, start, and
          end positions. Note, the BED format uses 0-based indexing, so make
          sure to number the first base as 0.
        - Sort the BED file by chromosome position (sort -k1,1 -k2,2n cpg.bed).
        - Gzip your sorted BED file and you now have your cpg.bed.gz QC file.

   2. windows100bp.gc_content.top10p.bed.gz
        - Group genome reference FASTA into 100 bp windows.
        - Calculate GC-content fraction for each window.
        - Create BED file with window positions (chromosome, start, and end) and
          GC-content fraction as a fourth column.
        - Sort BED file by GC-content fraction (sort -k4,4n gc_content.bed).
        - Find top 10% of GC-content windows.
        - Copy the four columns (chromosome, start, end, and GC-content
          fraction) for these windows into windows100bp.gc_content.top10p.bed.
        - Sort by position (sort -k1,1 -k2,2n) and gzip
          windows100bp.gc_content.top10p.bed to create your top 10% GC-content
          QC file.

   3. windows100bp.gc_content.bot10p.bed.gz
        - Follow Steps 1-4 for creating windows100bp.gc_content.top10p.bed.gz QC
          file.
        - Instead of finding the top 10% of GC-content windows, find the bottom
          10% of NON-zero GC-content windows.
        - Copy the four columns (chromosome, start, end, and GC-content
          fraction) for these windows into windows100bp.gc_content.bot10p.bed.
        - Sort by position (sort -k1,1 -k2,2n) and gzip
          windows100bp.gc_content.bot10p.bed to create your bottom 10%
          GC-content QC file.

Arguments:

    - Required:
        ref     - reference genome in FASTA format
        outdir  - output directory to save QC files to
    - Optional
        include - include all chromosomes/contigs
        verbose - print additional messages
        help    - print help/usage message

Dependencies:

    - command line: head, tail, sort, gzip
    - memory approximitating the size of your genome (whole reference is loaded)
=cut
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

# Usage message
my $usage = <<EOF;

Usage: perl build_biscuit_QC_assets.pl [-h,--help] [-v,--verbose] [-i,--include] <-o,--outdir directory> <-r,--ref reference>

Required inputs:
    -o,--outdir    : directory to put the 3 output biscuit QC asset files
    -r,--ref       : reference input genome in FASTA format
Optional flags:
    -i,--include   : include all chromosomes/contigs. Default is to restrict to chr[0-9]*,chrX,chrY,chrM
    -v,--verbose   : print additional messages
    -h,--help      : print usage

EOF

# Command line arguments
my($help,$ref,$outdir,$include,$verbose);
#=======================================================================
GetOptions(
    'o|outdir=s' => \$outdir, # string
    'r|ref=s' => \$ref,       # string
    'i|include' => \$include, # flag to include all chromosomes/contigs
    'v|verbose' => \$verbose, # flag for additional print statements
    'h|help' => \$help        # flag for help
);

# Check optional arguments
if (defined($help))    { print $usage; exit; }
if (defined($include)) { $include=1 } else { $include=0 }
if (defined($verbose)) { $verbose=1 } else { $verbose=0 }

# Check required arguments
if (!defined($outdir)) { die print "HELP: Define --outdir\n",$usage; }
if (!defined($ref))    { die print "HELP: Define --ref\n",$usage; }

# Setup directory for output files
if($outdir !~ /\/$/){
    $outdir = $outdir . '/';
    if ($verbose) { print STDOUT "$outdir changed (/ added to end)!\n"; }
}

if (! -e $outdir and ! -d $outdir) {
    print STDOUT "Creating output directory $outdir\n";
    mkdir $outdir
}

# Size of windows for calculating GC-content fraction
my $windowSize=100;

# Define outfiles
open(my $sw, '>', $outdir . "gc_content.bed") || die print "Cannot open gc_content.bed\n";

my $cpgBedName = "cpg.bed";
open(my $cpgBed, '>', "$outdir$cpgBedName") || die print "Cannot open $cpgBedName\n";

# Step 1. Load the reference genome
print STDOUT "Loading the reference genome\n";
my $in;
if ($ref =~ /.gz$/) {
    open($in, "gunzip -c $ref |") || die "can't open pipe to $ref\n";
} else {
    open($in, '<', $ref) || die "can't open $ref\n";
}
my $chr;
my %seq=();
while(<$in>) {
    chomp;
    if ($_ =~ /^>/ ) {
        #header
        $chr = $_;
        $chr =~ s/^>//;
        $chr =~ s/\s+.*//;
        if ($verbose) { print STDOUT "Found $chr\n"; }
    } else {
        # sequence
        $seq{$chr} .= $_;
    }
}

# Step 2. Find all CpGs
print STDOUT "Reference genome loaded, finding all CpGs\n";
my $char = 'CG';

foreach my $chr (sort keys %seq){
    # Check if we want to include all chromosomes/contigs
    if ($chr =~ /^chr\d{1,2}$/ or $chr=~/^chrM$/ or $chr=~/^chrY$/ or $chr=~/^chrX$/) {
        # Found a standard chromosome/contig - do nothing
        if ($verbose) { print STDOUT "Found a standard chr/contig: $chr\n"; }
    } else {
        # Found a nonstandard chromosome/contig
        if ($include) {
            # do nothing
            if ($verbose) { print STDOUT "Found a non-standard chr/contig: $chr - including\n"; }
        } else {
            # skip
            if ($verbose) { print STDOUT "Found a non-standard chr/contig: $chr - skipping\n"; }
            next;
        }
    }

    # Find the CpGs...
    my $offset = 0;          # Current offset/start position of CpG
    my $nCpGs=0;             # Number of CpGs
    my $str = uc $seq{$chr}; # Uppercase sequence for easier comparison

    my $result = index($str, $char, $offset);
    while ($result != -1) { # now get CpGs & print to bed
        #~ print "Found $char at $result\n";
        my $end = $result+2;
        print $cpgBed "$chr\t$result\t$end\n";
        $nCpGs++;
        $offset = $result + 1;
        $result = index($str, $char, $offset);
    }
    if ($verbose) { print STDOUT "\tchr/contig length = ", length($str), "    # CpGs = $nCpGs\n"; }
}
close($cpgBed);

# Step 3. Do the sliding window
print STDOUT "Sliding window analysis GC-content\n";
my $nWindows=0;

foreach my $chr (sort keys %seq){
    # Check if we want to include all chromosomes/contigs
    if ($chr =~ /^chr\d{1,2}$/ or $chr=~/^chrM$/ or $chr=~/^chrY$/ or $chr=~/^chrX$/) {
        # Found a standard chromosome/contig - do nothing
        if ($verbose) { print STDOUT "Found a standard chr/contig: $chr\n"; }
    } else {
        # Found a nonstandard chromosome/contig
        if ($include) {
            # do nothing
            if ($verbose) { print STDOUT "Found a non-standard chr/contig: $chr - including\n"; }
        } else {
            # skip
            if ($verbose) { print STDOUT "Found a non-standard chr/contig: $chr - skipping\n"; }
            next;
        }
    }

    # Move through each chromosome in windowSize bp chunks, calculating
    # GC-content for each chunk
    for (my $i=0; $i<=length($seq{$chr}); $i+=$windowSize) {
        my $substr = substr($seq{$chr}, $i, $windowSize); # chunk to search
        my @gcMatches = $substr =~ /[cg]/gi;              # number of Gs and Cs
        my @nMatches = $substr =~ /n/gi;                  # number of Ns

        my $gcFrac = sprintf('%.2f', (scalar @gcMatches/$windowSize));
        #my $gcFrac = sprintf('%f', (scalar @gcMatches/$windowSize));
        #if ($gcFrac == "0.0") {
        #    $gcFrac = "0"
        #} else {
        #    if ($gcFrac == "1.0") {
        #        $gcFrac = "1"
        #    } else {
        #        $gcFrac =~ s/0+$//; # remove trailing zeroes
        #    }
        #}
        #~ print STDOUT "Found a substr of length ", length($substr), " with GC-content: $gcFrac\n";
        #~ print STDOUT "\t$substr\n";

        if ((scalar @nMatches == 0) and (length($substr) == $windowSize)) {
            $nWindows++;
            my $end = $i+$windowSize;
            print $sw "$chr\t$i\t$end\t$gcFrac\n";
        }
    }
}

# Step 4. Command lines to sort & get the top/bottom 10% files...
print STDOUT "Command line functions to sort and subset files\n";

my $tenPerc = sprintf('%.0f',0.1*$nWindows);
if ($verbose) { print STDOUT "10% of $nWindows 100bp CpG windows is $tenPerc\n"; }

# Get top/bottom 10% GC-content windows, sort, and compress
system("LC_ALL=C sort -k4,4n $outdir/gc_content.bed > $outdir/gc_content.sorted.bed");
#system("LC_ALL=C sort -k4,4 -k1,1 -k2,2n $outdir/gc_content.bed > $outdir/gc_content.sorted.bed");
system("head -n $tenPerc $outdir/gc_content.sorted.bed | LC_ALL=C sort -k1,1 -k2,2n | gzip -c > $outdir/windows100bp.gc_content.bot10p.bed.gz");
system("tail -n $tenPerc $outdir/gc_content.sorted.bed | LC_ALL=C sort -k1,1 -k2,2n | gzip -c > $outdir/windows100bp.gc_content.top10p.bed.gz");

# Compress CpG output file
system("gzip $outdir/cpg.bed");

# Remove intermediate files
system("rm $outdir/gc_content.bed");
system("rm $outdir/gc_content.sorted.bed");

print "\nFinished $0\n";
#=======================================================================
#( Subroutines                  )
# ------------------------------------ 
#  o
#   o   \_\_    _/_/
#    o      \__/
#           (oo)\_______
#           (__)\       )\/\
#               ||----w |
#               ||     ||
