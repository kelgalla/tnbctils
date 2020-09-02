#/usr/bin/env perl

#---------------------------------------------------------------------
# FILE     : runCurl.pl
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-08-16
# COMMENTS : run curl queries to TCGA API
#---------------------------------------------------------------------

#curl 'https://api.gdc.cancer.gov/legacy/files/bfd9f7ad-f6bd-42dd-9687-dca1a605d35f?pretty=true'

#curl 'https://api.gdc.cancer.gov/legacy/files/ffe640d4-b79e-44d1-ab76-481bd2620561?pretty=true'

use URI::Escape;
use JSON;

use FindBin qw($Bin);
use lib "$Bin/..";
use ANALYSIS::ANALYSISShare qw(uuidToBarcode uuidToData uuidSampleType uuidToBarcodeFull);

$| = 1;

#my $uuidList = ["ffe640d4-b79e-44d1-ab76-481bd2620561", "bfd9f7ad-f6bd-42dd-9687-dca1a605d35f"];

my $uuidList = ["203a411d-21eb-4ff2-aabf-367303482f9b"];

#my $uuidList = ["5dd56b1e-be8e-46a2-9bd2-c7da4bb7e30d", "6c2aad88-4a1e-444a-a4f3-99250640170d"];

#my $res = &uuidToBarcode($uuidList);

#foreach my $uuid (keys %{$res}){
#    my $barcode = $res->{$uuid};
#    print($uuid . "\t" . $barcode . "\n");
#}

my $res = &uuidToBarcode($uuidList);

foreach my $uuid (keys %{$res}){
    my $type = $res->{$uuid};
    print($uuid . "\t" . $type . "\n");
}


#my @uuidList2 = map {'"' . $_ . '"'} @{$uuidList};
#my $uuidContent = join(",", @{$uuidList});

#my $json = {
#    "filters" => {
#	"op" => "in",
#	"content" => {
#	    "field" => "file_id",
#	    "value" => $uuidList
#	}},
#    "format" => "JSON",
#    "fields" => "cases.submitter_id",
#    "pretty" => "true"
#};

#$json = encode_json($json);
#print($json . "\n\n");

#my $json = '{"filters":{"op":"in","content":{"field":"file_id","value":[' . $uuidContent . ']}},';
#$json .= '"format":"JSON",';
#$json .= '"fields":"cases.submitter_id",';
#$json .= '"pretty":"true"}';
#print($json . "\n\n");

#$json = uri_escape($json);
#print($json . "\n\n");
#my $url = "https://api.gdc.cancer.gov/legacy/files?filters=" . $json . "&fields=cases.submitter_id&pretty=true";
#print($url . "\n\n");

#my $cmd = "curl --request POST --header \"Content-Type: application/json\" --data '$json' 'https://api.gdc.cancer.gov/legacy/files'";
#print($cmd . "\n\n");
#my $resjson = qx/$cmd/;
#print($resjson . "\n\n");

#my $res = decode_json($resjson);
#my $hits = $res->{"data"}->{"hits"}; # array ref

#foreach my $hit (@{$hits}){
#    my $id = $hit->{"id"};
#    my $barcode = $hit->{"cases"}->[0]->{"submitter_id"};
#    print($id . "\t" . $barcode . "\n");
#}



#my $cmd = `curl 'https://api.gdc.cancer.gov/legacy/files?filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22file_id%22%2C%22value%22%3A%5B%22bfd9f7ad-f6bd-42dd-9687-dca1a605d35f%22%2C%22ffe640d4-b79e-44d1-ab76-481bd2620561%22%5D%7D%7D&pretty=true'`;

#print($cmd . "\n");
