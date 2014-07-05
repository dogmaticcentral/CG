CREATE TABLE  aa_sequence (
	 ensembl_id  VARCHAR (20), 
	 peptide    BLOB
) ENGINE=MyISAM;

create index ensmbl_idx   on aa_sequence (ensembl_id);
