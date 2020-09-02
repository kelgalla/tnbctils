#/usr/bin/env perl
#---------------------------------------------------------------------
# FILE     : barcodeFromMAF.pl
# AUTHOR   : Kelly E. Craven
# CREATED  : 2018-02-18
# COMMENTS : isolate the unique barcodes from the maf file and compare 
#            to the tnbc cases to see if all cases have mutation info
#---------------------------------------------------------------------

use strict;
use warnings;

#use FindBin qw($Bin);
#use lib "$Bin/..";
#use ANALYSIS::ANALYSISShare qw(getBiotab);

# maf file location
my $mafFile = "/N/u/kelgalla/Karst/brca/brcatcga/dataReorg/maf/TCGA.BRCA.muse.d565f30b-d3bb-4739-b26c-be8d43c596d4.DR-6.0.protected.maf";

# tnbc file location
my $tnbcFile = "/N/u/kelgalla/Karst/brca/brcatcga/analysis/R/133tnbcCases.txt";

# read the unique maf barcodes into this hash
my %barcode;
open(MAF, "<$mafFile") or die("Couldn't open '$mafFile' for reading: ($!)\n");

my $header;
while ($header = <MAF>){
    if ($header =~ m/^\#/){
	next;
    } else {
	last;
    }
}
chomp($header);
my @header = split("\t", $header);
my %indices = map { $header[$_], $_} (0..$#header);

while (my $line = <MAF>){
    chomp($line);
    my @values = split("\t", $line);

    my $fullbcode = $values[$indices{"Tumor_Sample_Barcode"}];
    if ($fullbcode =~ m/(TCGA-[a-zA-Z0-9]{2}-[a-zA-Z0-9]{4})/) {
	my $bcode = $1;
	$barcode{$bcode} = 1;
    } else {
	die("Improper barcode\n");
    }
}

close(MAF);

#foreach my $key (keys %barcode){
#    print($key . "\n");
#}

# read the tnbc cases into this hash
my %tnbc;
open(TNBC, "<$tnbcFile") or die("Couldn't open '$tnbcFile' for reading: ($!)\n");

while (my $line = <TNBC>){
    chomp($line);
    $tnbc{$line} = "no mutations";
}

close(TNBC);

# mark if the tnbc cases are found in the maf file
foreach my $key (keys %barcode){
    if ($tnbc{$key}){
	$tnbc{$key} = "mutation found";
    }
}

foreach my $key (keys %tnbc){
    print($key . "\t" . $tnbc{$key} . "\n");
}
