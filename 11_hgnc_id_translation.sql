CREATE TABLE `id_translation` (
       `approved_symbol`  varchar(40)   DEFAULT NULL,
       `approved_name`    blob  DEFAULT NULL,
       `previous_symbols` blob  DEFAULT NULL,
`previous_names`   blob  DEFAULT NULL,
`synonyms`         blob  DEFAULT NULL,
`entrez_gene_id`   varchar(10)   DEFAULT NULL,
`ensembl_gene_id_hgnc`  varchar(40)   DEFAULT NULL,
`mouse_genome_database_id` varchar(40)   DEFAULT NULL,
`pubmed_ids`      blob  DEFAULT NULL,
`refseq_ids`      blob  DEFAULT NULL,
`uniprot_Ids`     blob  DEFAULT NULL,
`ensembl_gene_id`  varchar(40)   DEFAULT NULL,
  PRIMARY KEY (`approved_symbol`),
  KEY `approved_symbol_idx` (`approved_symbol`),
  KEY `entrez_idx` (`entrez_gene_id`)
) ENGINE=MyISAM;
