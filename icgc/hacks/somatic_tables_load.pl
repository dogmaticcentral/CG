#!/usr/bin/perl
#
# This source code is part of icgc, an ICGC processing pipeline.
#
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
#
# Contact: ivana.mihalek@gmail.com
#

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
    `mysql icgc < $table_name.sql`;
}


0;