#/usr/bin/env perl
#---------------------------------------------------------------------
# FILE     : isolateTNBC.pl
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-05-21
# COMMENTS : read the TNBC cases from the tnbc.txt file, and copy the 
#            relevant files to a separate directory
#---------------------------------------------------------------------

use strict;
use warnings;

use File::Copy;

#use FindBin qw($Bin);
#use lib "$Bin/..";
#use ANALYSIS::ANALYSISShare qw(uuidToBarcode);

# tnbc file location
my $tnbcFile = "/N/u/kelgalla/Karst/brca/analysis/R/tnbc.txt";

# normal
my $dataReorgDir = "/N/u/kelgalla/Karst/brca/tnbc";

# htseq
#my $dataTypeDir = "htseq";

# fpkm
#my $dataTypeDir = "fpkm";

# fpkm-uq
#my $dataTypeDir = "fpkm-uq";

# bcrXML
#my $dataTypeDir = "bcrXML";

# biospecimen
#my $dataTypeDir = "biospecimen";

# legacy
my $dataLegReorgDir = "/N/u/kelgalla/Karst/brca/dataLegacyReorg";

# clinSupp
#my $dataTypeDir = "clinSupp";
#$dataReorgDir = $dataLegReorgDir;

# pdf
#my $dataTypeDir = "pdf";
#$dataReorgDir = $dataLegReorgDir;

# svs
my $dataTypeDir = "svs";
$dataReorgDir = $dataLegReorgDir;

# tnbc
my $tnbcDir = "/N/u/kelgalla/Karst/brca/tnbc";

# read the barcodes from the file
my %barcode;
open(TNBC, "<$tnbcFile") or die("Couldn't open '$tnbcFile' for reading: ($!)\n");

my $header = <TNBC>;
chomp($header);
my @header = split("\t", $header);
my %indices = map { $header[$_], $_} (0..$#header);

while (my $line = <TNBC>){
    chomp($line);
    my @values = split("\t", $line);

    my $bcode = $values[$indices{"bcr_patient_barcode"}];
    $barcode{$bcode} = [];
}
close(TNBC);

# loop through files and add them to hash indexed by barcode
opendir(REORG, "$dataReorgDir/$dataTypeDir") or die("opendir() failed: $!");
while(defined(my $file = readdir(REORG))){
    if($file eq "." || $file eq ".."){next;}
    if ($file =~ m/(TCGA-[a-zA-Z0-9]{2}-[a-zA-Z0-9]{4})/){
	my $bcode = $1;
	if (defined($barcode{$bcode})) {
	    push(@{$barcode{$bcode}}, $file);
	}
    }
}
closedir(REORG);

# for each barcode, copy file(s) to new dir
my $path;
my $newPath;
foreach my $bcode (keys %barcode){
    if (@{$barcode{$bcode}}){
	print($bcode . ": " . scalar(@{$barcode{$bcode}}) . " files processed\n");
	foreach my $file (@{$barcode{$bcode}}){
	    $path = "$dataReorgDir/$dataTypeDir/$file";
	    $newPath = "$tnbcDir/$dataTypeDir/$file";
	    copy($path, $newPath) or die "Copy failed: $!";
	}
    } else {
	print($bcode . ": not found\n");
    }
}

