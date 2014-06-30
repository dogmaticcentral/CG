#!/usr/bin/perl

@maf_files = split '\n', `find . -name "*maf"`;

@common_fields = ();

foreach $maffile (@maf_files) {

    print "$maffile\n";
    open (IF, "<$maffile");
    while (<IF>) {
	next if ($_ =~ /^#/);
	@header_fields = map {lc} split;
	
	if (!@common_fields) {
	   @common_fields = @header_fields;
	   push @common_fields, 'aa_change';
	}

	@new_common_fields = ();
	foreach $hf (@header_fields) {
	    foreach $aa_change_name ('amino_acid_change_wu',
				     'aachange', 'amino_acid_change',
				     'protein_change'  ) {
		if ($hf eq $aa_change_name) {
		    print ">>>  $hf  $aa_change_name\n";
		    $hf = 'aa_change';

		    last;
		}
	    }
	    if ( grep {$hf eq $_} @common_fields) {
		push  @new_common_fields, $hf;
	    }
	}
	print "\n";
	@common_fields = @new_common_fields;
	

	last;
    }
    close IF;
}

#print join "\n", @common_fields;
#print "\n";
