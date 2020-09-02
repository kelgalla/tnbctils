#/usr/bin/env perl
#---------------------------------------------------------------------
# FILE     : condenseSurv.pl
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-08-04
# COMMENTS : read in a tab-delimited file with TCGA info, including 
#            survival information, and summarize it for analysis by R
#---------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Std;
use List::Util qw(max min);

use FindBin qw($Bin);
use lib "$Bin/..";
use ANALYSIS::ANALYSISShare qw();

$| = 1;

our ($opt_i, $opt_o);
my $usage = "Read TCGA data file and add a status and days column to summarize the survival data\n";
$usage .= "perl condenseSurv.pl -i tcgaFile\n";
$usage .= "-i file    input tcga file with survival data after running extractInfo.pl\n";
$usage .= "-o file    output tcga file with a status and days column of summarized survival data\n\n";

getopts("i:o:") or die($usage);

if (!$opt_i or !-e $opt_i){
    die("condenseSurv.pl: You must have a valid input file.\n$usage");
}
if (!$opt_o){
    die("condenseSurv.pl: You must specify an output file.\n$usage");
}

open(IN, "<$opt_i") or die("Counldn't open '$opt_i' for reading: ($!)\n");
open(OUT, ">$opt_o") or die("Couldn't open '$opt_o' for writing: ($!)\n");

my $header = <IN>;
my @header = split("\t", $header);
chomp(@header);
push(@header, ("os_status", "os_days", "dfs_status", "dfs_days"));
print(OUT join("\t", @header), "\n");
my %indices = map { $header[$_], $_} (0..$#header);

while (my $line = <IN>){
    my @surv; # (censored (0) or dead(1), days)
    my @dfs; # (censored (0) or dead/nte(1), days)

    my @values = split("\t", $line);
    chomp(@values);

    my $bcode = $values[$indices{"bcr_patient_barcode"}];

    # initial data
    my $vital_status = $values[$indices{"vital_status"}];
    my @vital_status = split(/ \| /, $vital_status);
    my $days_to_death = $values[$indices{"days_to_death"}];
    my @days_to_death = split(/ \| /, $days_to_death);
    my $days_to_last_followup = $values[$indices{"days_to_last_followup"}];
    my @days_to_last_followup = split(/ \| /, $days_to_last_followup);

    # v1.5 data
    my $one_five_vital_status = $values[$indices{"v1.5_vital_status"}];
    my @one_five_vital_status = split(/ \| /, $one_five_vital_status);
    my $one_five_days_to_death = $values[$indices{"v1.5_days_to_death"}];
    my @one_five_days_to_death = split(/ \| /, $one_five_days_to_death);
    my $one_five_days_to_last_followup = $values[$indices{"v1.5_days_to_last_followup"}];
    my @one_five_days_to_last_followup = split(/ \| /, $one_five_days_to_last_followup);

    # v2.1 data
    my $two_one_vital_status = $values[$indices{"v2.1_vital_status"}];
    my @two_one_vital_status = split(/ \| /, $two_one_vital_status);
    my $two_one_days_to_death = $values[$indices{"v2.1_days_to_death"}];
    my @two_one_days_to_death = split(/ \| /, $two_one_days_to_death);
    my $two_one_days_to_last_followup = $values[$indices{"v2.1_days_to_last_followup"}];
    my @two_one_days_to_last_followup = split(/ \| /, $two_one_days_to_last_followup);

    # v4.0 data
    my $four_zero_vital_status = $values[$indices{"v4.0_vital_status"}];
    my @four_zero_vital_status = split(/ \| /, $four_zero_vital_status);
    my $four_zero_days_to_death = $values[$indices{"v4.0_days_to_death"}];
    my @four_zero_days_to_death = split(/ \| /, $four_zero_days_to_death);
    my $four_zero_days_to_last_followup = $values[$indices{"v4.0_days_to_last_followup"}];
    my @four_zero_days_to_last_followup = split(/ \| /, $four_zero_days_to_last_followup);


    # NTE section
    my $new_tumor_event_after_initial_treatment =
	$values[$indices{"new_tumor_event_after_initial_treatment"}];
    
    my $days_to_new_tumor_event_after_initial_treatment = 
	$values[$indices{"days_to_new_tumor_event_after_initial_treatment"}];
    my @days_to_new_tumor_event_after_initial_treatment = 
	split(/ \| /, $days_to_new_tumor_event_after_initial_treatment);

    my $days_to_new_tumor_event_additional_surgery_procedure = 
	$values[$indices{"days_to_new_tumor_event_additional_surgery_procedure"}];
    my @days_to_new_tumor_event_additional_surgery_procedure = 
	split(/ \| /, $days_to_new_tumor_event_additional_surgery_procedure);

    # combine the 2 fields into one
    my @days_to_new_tumor_event1;
    foreach my $i (0..$#days_to_new_tumor_event_after_initial_treatment){
	if ($days_to_new_tumor_event_after_initial_treatment[$i] =~ /\d+/){
	    push(@days_to_new_tumor_event1, $days_to_new_tumor_event_after_initial_treatment[$i]);
	} elsif ($days_to_new_tumor_event_additional_surgery_procedure[$i] =~ /\d+/) {
	    push(@days_to_new_tumor_event1, $days_to_new_tumor_event_additional_surgery_procedure[$i]);
	}
    }

    my $new_neoplasm_event_type =
	$values[$indices{"new_neoplasm_event_type"}];
    my @new_neoplasm_event_type = 
	split(/ \| /, $new_neoplasm_event_type);

    my @days_to_new_tumor_event;
    for (0..$#new_neoplasm_event_type){
	if ($new_neoplasm_event_type[$_] ne "New Primary Tumor"){
	    push(@days_to_new_tumor_event, 
		 $days_to_new_tumor_event1[$_]);
	}
    }

    # v1.5
    my $one_five_days_to_new_tumor_event_after_initial_treatment = 
	$values[$indices{"v1.5_days_to_new_tumor_event_after_initial_treatment"}];
    my @one_five_days_to_new_tumor_event = 
	split(/ \| /, $one_five_days_to_new_tumor_event_after_initial_treatment);

    # v2.1
    my $two_one_new_tumor_event_after_initial_treatment = 
	$values[$indices{"v2.1_new_tumor_event_after_initial_treatment"}];
    my @two_one_new_tumor_event_after_initial_treatment = 
	split(/ \| /, $two_one_new_tumor_event_after_initial_treatment);

    my $two_one_days_to_new_tumor_event_after_initial_treatment = 
	$values[$indices{"v2.1_days_to_new_tumor_event_after_initial_treatment"}];
    my @two_one_days_to_new_tumor_event_after_initial_treatment = 
	split(/ \| /, $two_one_days_to_new_tumor_event_after_initial_treatment);

    my @two_one_days_to_new_tumor_event;
    if (my $nteNum = grep {$_ eq "YES" } @two_one_new_tumor_event_after_initial_treatment){
	if (@two_one_days_to_new_tumor_event_after_initial_treatment != $nteNum){
	    foreach my $i (0..$#two_one_new_tumor_event_after_initial_treatment){
		if ($two_one_new_tumor_event_after_initial_treatment[$i] eq "YES"){
		    if ($two_one_vital_status[$i] eq "Alive"){
			push(@two_one_days_to_new_tumor_event, $two_one_days_to_last_followup[$i]);
		    } elsif ($two_one_vital_status[$i] eq "Dead"){
			push(@two_one_days_to_new_tumor_event, $two_one_days_to_death[$i]);
		    }
		}
	    }
	}
    }

    my $two_one_new_neoplasm_event_type =
	$values[$indices{"v2.1_new_neoplasm_event_type"}];
    my @two_one_new_neoplasm_event_type = 
	split(/ \| /, $two_one_new_neoplasm_event_type);

    for (0..$#two_one_new_neoplasm_event_type){
	if ($two_one_new_neoplasm_event_type[$_] ne "New Primary Tumor"){
	    push(@two_one_days_to_new_tumor_event, 
		 $two_one_days_to_new_tumor_event_after_initial_treatment[$_]);
	}
    }
    
    # v4.0
    my $four_zero_new_tumor_event_after_initial_treatment = 
	$values[$indices{"v4.0_new_tumor_event_after_initial_treatment"}];
    my @four_zero_new_tumor_event_after_initial_treatment = 
	split(/ \| /, $four_zero_new_tumor_event_after_initial_treatment);

    my $four_zero_days_to_new_tumor_event_after_initial_treatment = 
	$values[$indices{"v4.0_days_to_new_tumor_event_after_initial_treatment"}];
    my @four_zero_days_to_new_tumor_event_after_initial_treatment = 
	split(/ \| /, $four_zero_days_to_new_tumor_event_after_initial_treatment);

    my $four_zero_new_neoplasm_event_type =                       
	$values[$indices{"v4.0_new_neoplasm_event_type"}];
    my @four_zero_new_neoplasm_event_type = 
	split(/ \| /, $four_zero_new_neoplasm_event_type);

    my @four_zero_days_to_new_tumor_event;
    for (0..$#four_zero_new_neoplasm_event_type){
	if ($four_zero_new_neoplasm_event_type[$_] ne "New Primary Tumor"){
	    push(@four_zero_days_to_new_tumor_event, 
		 $four_zero_days_to_new_tumor_event_after_initial_treatment[$_]);
	}
    }

    my $nteDFS;
    my $nteSURV;
    if ((@days_to_new_tumor_event >= 1) || (@one_five_days_to_new_tumor_event >=1) || 
	(@two_one_days_to_new_tumor_event >=1) || (@four_zero_days_to_new_tumor_event >=1)){
	my @newarray = (@days_to_new_tumor_event, @one_five_days_to_new_tumor_event,
			@two_one_days_to_new_tumor_event, @four_zero_days_to_new_tumor_event);
	@newarray = grep { $_ ne "" } @newarray;
	$nteDFS = min(@newarray);
	$nteSURV = max(@newarray);
    }
    if (defined($nteDFS)){
	@dfs = (1, $nteDFS);
    }
    if (defined($nteSURV)){
	@surv = (0, $nteSURV);
    }


    # initial data
    if (@vital_status == 0){
	# do nothing
    } elsif (@vital_status >= 1){
	foreach my $i (0..$#vital_status){
	    if ($vital_status[$i] eq "Alive"){
		if (@surv == 0){ # not set yet
		    @surv = (0, $days_to_last_followup[$i]);
		} elsif (@surv && ($surv[0] == 0)){ # already set, determine which number is higher
		    my $survDays = $surv[1];
		    my $days = $days_to_last_followup[$i];
		    @surv = (0, max($survDays, $days));
		}
		if (@dfs == 0){ # not set yet
		    @dfs = (0, $days_to_last_followup[$i]);
		} elsif (@dfs && ($dfs[0] == 0)){
		    my $dfsDays = $dfs[1];
		    my $days = $days_to_last_followup[$i];
		    @dfs = (0, max($dfsDays, $days));
		}
	    } else { # dead
		if (@surv == 0){ # not set yet
		    @surv = (1, $days_to_death[$i]);
		} elsif (@surv && ($surv[0] == 0)) { # this will trump it
		    @surv = (1, $days_to_death[$i]);
		} elsif (@surv && ($surv[0] == 1)){
		    my $survDays = $surv[1];
		    my $days = $days_to_death[$i];
		    @surv = (1, max($survDays, $days));
		}
		if (@dfs == 0){ # not set yet
		    @dfs = (1, $days_to_death[$i]);
		} elsif (@dfs && ($dfs[0] == 0)){ # this will trump it
		    @dfs = (1, $days_to_death[$i]);
		} elsif (@dfs && ($dfs[0] == 1)){
		    # do nothing
		}
	    }
	}
    }


    # v1.5 data
    if (@one_five_vital_status == 0){
	# do nothing
    } elsif (@one_five_vital_status >= 1){
	foreach my $i (0..$#one_five_vital_status){
	    if ($one_five_vital_status[$i] eq "Alive"){
		if (@surv == 0){ # not set yet
		    @surv = (0, $one_five_days_to_last_followup[$i]);
		} elsif (@surv && ($surv[0] == 0)){ # already set, determine which number is higher
		    my $survDays = $surv[1];
		    my $days = $one_five_days_to_last_followup[$i];
		    @surv = (0, max($survDays, $days));
		}
		if (@dfs == 0){ # not set yet
		    @dfs = (0, $one_five_days_to_last_followup[$i]);
		} elsif (@dfs && ($dfs[0] == 0)){
		    my $dfsDays = $dfs[1];
		    my $days = $one_five_days_to_last_followup[$i];
		    @dfs = (0, max($dfsDays, $days));
		}
	    } else { # dead
		if (@surv == 0){ # not set yet
		    @surv = (1, $one_five_days_to_death[$i]);
		} elsif (@surv && ($surv[0] == 0)) { # this will trump it
		    @surv = (1, $one_five_days_to_death[$i]);
		} elsif (@surv && ($surv[0] == 1)){
		    my $survDays = $surv[1];
		    my $days = $one_five_days_to_death[$i];
		    @surv = (1, max($survDays, $days));
		}
		if (@dfs == 0){ # not set yet
		    @dfs = (1, $one_five_days_to_death[$i]);
		} elsif (@dfs && ($dfs[0] == 0)){ # this will trump it
		    @dfs = (1, $one_five_days_to_death[$i]);
		} elsif (@dfs && ($dfs[0] == 1)){
		    # do nothing
		}
	    }
	}
    }


    # v2.1 data
    if (@two_one_vital_status == 0){
	# do nothing
    } elsif (@two_one_vital_status >= 1){
	foreach my $i (0..$#two_one_vital_status){
	    if ($two_one_vital_status[$i] eq "Alive"){
		if (@surv == 0){ # not set yet
		    @surv = (0, $two_one_days_to_last_followup[$i]);
		} elsif (@surv && ($surv[0] == 0)){ # already set, determine which number is higher
		    my $survDays = $surv[1];
		    my $days = $two_one_days_to_last_followup[$i];
		    @surv = (0, max($survDays, $days));
		}
		if (@dfs == 0){ # not set yet
		    @dfs = (0, $two_one_days_to_last_followup[$i]);
		} elsif (@dfs && ($dfs[0] == 0)){
		    my $dfsDays = $dfs[1];
		    my $days = $two_one_days_to_last_followup[$i];
		    @dfs = (0, max($dfsDays, $days));
		}
	    } else { # dead
		if (@surv == 0){ # not set yet
		    @surv = (1, $two_one_days_to_death[$i]);
		} elsif (@surv && ($surv[0] == 0)) { # this will trump it
		    @surv = (1, $two_one_days_to_death[$i]);
		} elsif (@surv && ($surv[0] == 1)){
		    my $survDays = $surv[1];
		    my $days = $two_one_days_to_death[$i];
		    @surv = (1, max($survDays, $days));
		}
		if (@dfs == 0){ # not set yet
		    @dfs = (1, $two_one_days_to_death[$i]);
		} elsif (@dfs && ($dfs[0] == 0)){ # this will trump it
		    @dfs = (1, $two_one_days_to_death[$i]);
		} elsif (@dfs && ($dfs[0] == 1)){
		    # do nothing
		}
	    }
	}
    }


    # v4.0 data
    if (@four_zero_vital_status == 0){
	# do nothing
    } elsif (@four_zero_vital_status >= 1){
	foreach my $i (0..$#four_zero_vital_status){
	    if ($four_zero_vital_status[$i] eq "Alive"){
		if (@surv == 0){ # not set yet
		    @surv = (0, $four_zero_days_to_last_followup[$i]);
		} elsif (@surv && ($surv[0] == 0)){ # already set, determine which number is higher
		    my $survDays = $surv[1];
		    my $days = $four_zero_days_to_last_followup[$i];
		    @surv = (0, max($survDays, $days));
		}
		if (@dfs == 0){ # not set yet
		    @dfs = (0, $four_zero_days_to_last_followup[$i]);
		} elsif (@dfs && ($dfs[0] == 0)){
		    my $dfsDays = $dfs[1];
		    my $days = $four_zero_days_to_last_followup[$i];
		    @dfs = (0, max($dfsDays, $days));
		}
	    } else { # dead
		if (@surv == 0){ # not set yet
		    @surv = (1, $four_zero_days_to_death[$i]);
		} elsif (@surv && ($surv[0] == 0)) { # this will trump it
		    @surv = (1, $four_zero_days_to_death[$i]);
		} elsif (@surv && ($surv[0] == 1)){
		    my $survDays = $surv[1];
		    my $days = $four_zero_days_to_death[$i];
		    @surv = (1, max($survDays, $days));
		}
		if (@dfs == 0){ # not set yet
		    @dfs = (1, $four_zero_days_to_death[$i]);
		} elsif (@dfs && ($dfs[0] == 0)){ # this will trump it
		    @dfs = (1, $four_zero_days_to_death[$i]);
		} elsif (@dfs && ($dfs[0] == 1)){
		    # do nothing
		}
	    }
	}
    }

    push(@values, @surv, @dfs);
    print(OUT join("\t", @values), "\n");
}
close(OUT);
close(IN);
