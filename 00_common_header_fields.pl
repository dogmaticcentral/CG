#!/usr/bin/perl

$tcga_dir = "/Users/ivana/databases/TCGA/";
$cmd = 


@maf_files = split '\n', `find $tcga_dir -name "*maf"`;

@common_fields = ();

foreach $maffile (@maf_files) {

    #print "$maffile\n";
    open (IF, "<$maffile");
    while (<IF>) {
	next if ($_ =~ /^#/);
	@header_fields = map {lc} split;

	if (!@common_fields) {
	   @common_fields = @header_fields;
	   push @common_fields, 'aa_change';
	   push @common_fields, 'cdna_change';
	}

	@new_common_fields = ();
	$found = 0;
	$aa_name_used = 'none';
	foreach $hf (@header_fields) {
	    # give a common name to aa change field
	    foreach $aa_change_name ('amino_acid_change_wu',
				     'aachange', 'amino_acid_change',
				     'protein_change'  ) {
		if ($hf eq $aa_change_name) {
		    #print "###   $aa_change_name\n";
		    $aa_name_used = $hf;
		    $hf = 'aa_change';

		    last;
		}
	    }

	    foreach $cdna_change_name ('cdna_change', 'chromchange', 'c_position_wu', 'c_position' ) {
		if ($hf eq $cdna_change_name) {
		    #print ">>>   $cdna_change_name\n";
		    $hf = 'cdna_change';
		    $found = 1;
		    last;
		}
	    }
	    if ( grep {$hf eq $_} @common_fields) {
		push  @new_common_fields, $hf;
	    }
	}
	if (!$found) {
	    #print "$maffile      $aa_name_used\n";
	    @new_common_fields = @common_fields;
	}

	#print "\n";
	@common_fields = @new_common_fields;

	last;
    }
    close IF;
}

print join "\n", @common_fields;
print "\n";
