#!/usr/bin/perl

use warnings;
use strict;


my %pos = ();
while (<>) {

    my @fields = split '\t';
    $fields[9] eq "pathogenic" && next;
    my $mut = $fields[7];
    defined $fields[7]|| next;
    $mut =~ /(\D+)(\d+)(\D+)/;
    defined $1 || next;
    defined $2 || next;
    defined $3 || next;
    $3 =~ /[\*\?]/  && next;
    $1 =~ /[\*\?]/  && next;
    length($1)==1 &&  length($3)==1 || next;
    #print "$fields[9] \n"; # $mut   $1   $2  $3\n";
    $pos{$2} = $1;

}

for my $p (sort {$a <=> $b} keys %pos) {
    print "$pos{$p}    $p\n";
}