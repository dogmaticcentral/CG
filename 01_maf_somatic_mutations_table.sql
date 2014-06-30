  CREATE TABLE somatic_mutations (
  	 hugo_symbol VARCHAR (50) NOT NULL, 

	 aa_change BLOB, 

	 chromosome VARCHAR (20) NOT NULL, 
	 start_position INT  NOT NULL, 
	 end_position INT NOT NULL, 
	 strand VARCHAR (5) NOT NULL, 

	 variant_classification VARCHAR (50) NOT NULL, 
	 variant_type VARCHAR (20) NOT NULL, 

	 reference_allele BLOB NOT NULL, 
	 tumor_seq_allele1 BLOB NOT NULL, 
	 tumor_seq_allele2 BLOB NOT NULL, 

	 tumor_sample_barcode VARCHAR (50) NOT NULL, 
	 matched_norm_sample_barcode VARCHAR (50) NOT NULL, 

	 match_norm_seq_allele1 BLOB, 
	 match_norm_seq_allele2 BLOB, 
	 tumor_validation_allele1 BLOB, 
	 tumor_validation_allele2 BLOB, 
	 match_norm_validation_allele1 BLOB, 
	 match_norm_validation_allele2 BLOB, 

	 verification_status VARCHAR (20), 
	 validation_status VARCHAR (20) NOT NULL, 
	 mutation_status VARCHAR (50) NOT NULL

) ENGINE=MyISAM;


create index hugo_idx     on somatic_mutations (hugo_symbol);
create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position);

