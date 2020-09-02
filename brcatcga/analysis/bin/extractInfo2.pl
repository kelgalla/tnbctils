#/usr/bin/env perl
#---------------------------------------------------------------------
# FILE     : extractInfo2.pl
# AUTHOR   : Kelly E. Craven
# CREATED  : 2019-12-19
# COMMENTS : read the tnbc.tils.append.surv.idc.txt file, and append
#            additional information from the clinical bcrXML files
#---------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Std;

use FindBin qw($Bin);
use lib "$Bin/..";
use ANALYSIS::ANALYSISShare qw(getFields);

use File::Slurp;

$| = 1;

my $bcrXMLdir = "F://TNBC TILS/brcatcga/tnbc/bcrXML";
#my $bcrXMLdir = "/N/u/kelgalla/Karst/brca/brcatcga/tnbc/bcrXML";
#my $tnbctilsFile = "/N/u/kelgalla/Karst/brca/brcatcga/analysis/R/tnbc.tils.txt";
#my $tnbctilsappendFile = "/N/u/kelgalla/Karst/brca/brcatcga/analysis/R/tnbc.tils.append.txt";

our ($opt_i, $opt_o);
my $usage = "Read TCGA data file and add various columns of information including survival\n";
$usage .= "perl extractInfo.pl -i tcgaFile\n";
$usage .= "-i file    input tcga file after the tils information was added to the biotab file\n";
$usage .= "-o file    output tcga file with additional information including survival data\n\n";

getopts("i:o:") or die($usage);

if (!$opt_i or !-e $opt_i){
    die("extractInfo.pl: You must have a valid input file.\n$usage");
}
if (!$opt_o){
    die("extractInfo.pl: You must specify an output file.\n$usage");
}

# read the barcodes from the file
my %barcode;
open(TNBC, "<$opt_i") or die("Counldn't open '$opt_i' for reading: ($!)\n");

my $header = <TNBC>;
my @header = split("\t", $header);
chomp(@header);
my %indices = map { $header[$_], $_} (0..$#header);

while (my $line = <TNBC>){
    my @values = split("\t", $line);
    chomp(@values);
    
    my $bcode = $values[$indices{"bcr_patient_barcode"}];
    $barcode{$bcode} = [];
}
close(TNBC);

# loop through files and add them to hash indexed by barcode
opendir(REORG, "$bcrXMLdir") or die("opendir() failed: $!");
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

# loop through the tnbc tils files again, appending the relevant info
open(TNBC, "<$opt_i") or die("Counldn't open '$opt_i' for reading: ($!)\n");
open(TNBCAPP, ">$opt_o") or die("Couldn't open '$opt_o' for writing: ($!)\n");
my $header2 = <TNBC>;
my @header2 = split("\t", $header2);
chomp(@header2);
my @fields = ([q{*[local-name()="patient"]/*[local-name()="race_list"]} => ["race"]],
	      [q{*[local-name()="patient"]} => ["ethnicity", "menopause_status", "breast_cancer_surgery_margin_status", "axillary_lymph_node_stage_method_type", "axillary_lymph_node_stage_other_method_descriptive_text", "lymph_node_examined_count", "number_of_lymphnodes_positive_by_ihc", "number_of_lymphnodes_positive_by_he"]],
	      [q{*[local-name()="patient"]/*[local-name()="stage_event"]/*[local-name()="tnm_categories"]/*[local-name()="pathologic_categories"]} => ["pathologic_T", "pathologic_N", "pathologic_M"]],
	      [q{*[local-name()="patient"]/*[local-name()="first_nonlymph_node_metastasis_anatomic_sites"]} => ["metastatic_site_at_diagnosis", "metastatic_site_at_diagnosis_other"]],
	      [q{*[local-name()="patient"]} => ["distant_metastasis_present_ind2"]]);
# my @fields = ([q{*[local-name()="patient"]/*[local-name()="stage_event"]} => ["pathologic_stage"]],
# 	      [q{*[local-name()="patient"]} => ["age_at_initial_pathologic_diagnosis", "histological_type", "histological_type_other", "breast_carcinoma_surgical_procedure_name", "surgical_procedure_purpose_other_text", "margin_status", "breast_carcinoma_primary_surgical_procedure_name", "breast_neoplasm_other_surgical_procedure_descriptive_text"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="new_tumor_events"]} => ["new_tumor_event_after_initial_treatment"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="new_tumor_events"]/*[local-name()="new_tumor_event"]} => ["days_to_new_tumor_event_after_initial_treatment", "days_to_new_tumor_event_additional_surgery_procedure", "new_neoplasm_event_type"]],
# 	      [q{*[local-name()="patient"]} => ["days_to_death", "days_to_last_followup", "vital_status", "person_neoplasm_cancer_status"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="follow_ups"]/*[local-name()="follow_up"][@*[local-name()="version" and .="1.5"]]} => ["days_to_death", "days_to_last_followup", "vital_status", "person_neoplasm_cancer_status", "days_to_new_tumor_event_after_initial_treatment"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="follow_ups"]/*[local-name()="follow_up"][@*[local-name()="version" and .="2.1"]]} => ["days_to_death", "days_to_last_followup", "vital_status", "person_neoplasm_cancer_status", "new_tumor_event_after_initial_treatment", "days_to_new_tumor_event_after_initial_treatment", "new_neoplasm_event_type"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="follow_ups"]/*[local-name()="follow_up"][@*[local-name()="version" and .="4.0"]]} => ["days_to_death", "days_to_last_followup", "vital_status", "person_neoplasm_cancer_status"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="follow_ups"]/*[local-name()="follow_up"][@*[local-name()="version" and .="4.0"]]/*[local-name()="new_tumor_events"]} => ["new_tumor_event_after_initial_treatment"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="follow_ups"]/*[local-name()="follow_up"][@*[local-name()="version" and .="4.0"]]/*[local-name()="new_tumor_events"]/*[local-name()="new_tumor_event"]} => ["days_to_new_tumor_event_after_initial_treatment", "new_neoplasm_event_type"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="drugs"]/*[local-name()="drug"]} => ["days_to_drug_therapy_start", "days_to_drug_therapy_end", "drug_name"]],
# 	      [q{*[local-name()="patient"]/*[local-name()="radiations"]/*[local-name()="radiation"]} => ["days_to_radiation_therapy_start", "days_to_radiation_therapy_end", "radiation_dosage", "units"]],);

foreach my $tag (@fields){
    my $xpath = $tag->[0];
    my @list = @{$tag->[1]};
    if ($xpath =~ /version.*?(\d\.\d)/){
	@list = map { "v$1_" . $_ } @list;
    }
    push(@header2, @list);    
}
print(TNBCAPP join("\t", @header2), "\n");
my %indices2 = map { $header2[$_], $_} (0..$#header2);

while (my $line2 = <TNBC>){
    my @values2 = split("\t", $line2);
    chomp(@values2);
    my $bcode = $values2[$indices2{"bcr_patient_barcode"}];

    if (@{$barcode{$bcode}} == 1){
	my $file = $barcode{$bcode}->[0];
	my $xml = read_file("$bcrXMLdir/$file");
	foreach my $tag (@fields){
	    my $res = getFields($xml, $tag->[0], $tag->[1]);

	    my @fieldVals;
	    foreach my $field (@{$tag->[1]}){
		push(@fieldVals, $res->{$field});
	    }
	    push(@values2, @fieldVals);
	}
	print(TNBCAPP join("\t", @values2), "\n");
    } else {
	die("More or less than one bcrXML file associated with barcode $bcode: ($!)\n");
    }
}
close(TNBCAPP);
close(TNBC);
