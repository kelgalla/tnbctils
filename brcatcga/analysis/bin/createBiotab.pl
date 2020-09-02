#/usr/bin/env perl
#---------------------------------------------------------------------
# FILE     : isolateTNBC.pl
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-08-02
# COMMENTS : read the clinical bcrXML files and create a biotab.txt 
#            file in this directory
#---------------------------------------------------------------------

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/..";
use ANALYSIS::ANALYSISShare qw(getBiotab);

use File::Slurp;

my $bcrXMLdir = "/N/u/kelgalla/Karst/brca/brcatcga/dataReorg/bcrXML";
my $biotabFile = "/N/u/kelgalla/Karst/brca/brcatcga/analysis/bin/biotab.txt";

# find all clinical bcrXML files and put the name in the @xmlFile array
my @xmlFile;
opendir(BCRXML, "$bcrXMLdir") or die("opendir() failed: $!");
while(defined(my $file = readdir(BCRXML))){
    if($file eq "." || $file eq ".."){next;}
    push(@xmlFile, $file);
}
closedir(BCRXML);

open(BIOTAB, ">$biotabFile") or die("Couldn't open '$biotabFile' for writing: ($!)\n");

print(BIOTAB "bcr_patient_uuid\tbcr_patient_barcode" . "\t" . "breast_carcinoma_estrogen_receptor_status" . "\t" .
    "er_level_cell_percentage_category" . "\t" . "breast_carcinoma_progesterone_receptor_status" . "\t" .
    "progesterone_receptor_level_cell_percent_category" . "\t" . "lab_proc_her2_neu_immunohistochemistry_receptor_status" . "\t" .
    "her2_erbb_pos_finding_cell_percent_category" . "\t" . "her2_immunohistochemistry_level_result" . "\t" .
    "lab_procedure_her2_neu_in_situ_hybrid_outcome_type" . "\n");

# loop through the bcrXML files, isolate the biotab info, and write them to a file
foreach my $xmlFile (@xmlFile){

    my $xml = read_file("$bcrXMLdir/$xmlFile");
    my $res = getBiotab($xml);

    my $tbl .= $res->{"bcr_patient_uuid"} . "\t" . $res->{"bcr_patient_barcode"} . "\t" . $res->{"breast_carcinoma_estrogen_receptor_status"} . "\t" . 
	$res->{"er_level_cell_percentage_category"} . "\t" . $res->{"breast_carcinoma_progesterone_receptor_status"} . "\t" .
	$res->{"progesterone_receptor_level_cell_percent_category"} . "\t" . $res->{"lab_proc_her2_neu_immunohistochemistry_receptor_status"} . "\t" . 
	$res->{"her2_erbb_pos_finding_cell_percent_category"} . "\t" . $res->{"her2_immunohistochemistry_level_result"} . "\t" .
	$res->{"lab_procedure_her2_neu_in_situ_hybrid_outcome_type"} . "\n";

    print(BIOTAB $tbl);
}

close(BIOTAB);
