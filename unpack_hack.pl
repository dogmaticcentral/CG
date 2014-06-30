#!/usr/bin/perl 


@all = split "\n", `ls`;
$home = `pwd`; chomp $home;


@all_dirs = grep { -d $_ } @all;

foreach $dir  (@all_dirs) {
    chdir $home;
    chdir $dir;
    @tars =  split "\n", `ls *.tar`;
    foreach $tar (@tars) {
	`tar -xvf $tar && rm $tar`;
    }
    `rm -f */*.wig.*`;

}
