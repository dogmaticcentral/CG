#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

# grep -n "Table structure" icgc.May13.no_tcga.sql  > icgc.May13.no_tcga.table_bdries.txt
# find the number of lines in icgc.May13.no_tcga.sql somehow
# wc -l takes forever
# for example from the above icgc.May13.no_tcga.table_bdries.txt
# we know that 84461:-- Table structure for table `mutations_chrom_Y`
# and  tail -n 46 icgc.May13.no_tcga.sql  | head -n1
# ergo the number of lines is 84461+45=84506
# (this is crap. use as a last resort.)
my $number_of_lines = 84506;
my $infile  = "icgc.May13.no_tcga.table_bdries.txt";
my $dumpfile = "icgc.May13.no_tcga.sql";

my %linefrom = ();
my %lineto = ();
my @table_names;

open(INF,"<$infile") || die "CNo $infile: $!\n";
while(<INF>) {
    /somatic/ || next;
    chomp;
    my @aux = split ':';
    my $lineno = shift @aux;
    @aux = split;
    my $table_name = pop @aux; $table_name =~ s/\`//g;
    if ($table_name =~ /temp/) {
        $table_name =~ s/_temp//;
        $lineto{$table_name} = $lineno-1;
    } else {
        push @table_names, $table_name;
        $linefrom{$table_name} = $lineno;
    }

}
close INF;


foreach my $table_name (@table_names) {

    print "$table_name, $linefrom{$table_name}, $lineto{$table_name} \n";
    my $chunksize = $lineto{$table_name}  - $linefrom{$table_name} + 1;
    if ($linefrom{$table_name}<$number_of_lines/2 ) {
        `head -n $lineto{$table_name} $dumpfile | tail -n $chunksize > $table_name.sql`;
    } else {
        my $backward_count = $number_of_lines-$linefrom{$table_name}+1;
        `tail -n $backward_count $dumpfile | head -n $chunksize > $table_name.sql`
    }
}


0;