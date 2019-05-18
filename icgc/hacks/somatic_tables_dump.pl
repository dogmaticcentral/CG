#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

my $infile  = "icgc.May13.no_tcga.table_bdries.txt";
my @table_names;

open(INF,"<$infile") || die "CNo $infile: $!\n";
while(<INF>) {
    chomp;
     /somatic`$/ || next;
    my @aux = split;
    my $table_name = pop @aux; $table_name =~ s/\`//g;
    push @table_names, $table_name;
}
close INF;

foreach my $table_name (@table_names) {
    print "$table_name \n";
    `mysqldump icgc  $table_name > $table_name.sql`;
}


0;