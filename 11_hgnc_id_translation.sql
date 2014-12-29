/* this refers to the  hgnc_id_translation.txt table in ~/database */
/* it is supposed to go as on of the tables in basline database */
CREATE TABLE `hgnc_id_translation` (
       `hgnc_id` int NOT NULL,
       `approved_symbol`  varchar(40)   DEFAULT NULL,
       `approved_name`    blob  DEFAULT NULL,
 `locus_type`   blob  DEFAULT NULL,
      `previous_symbols` blob  DEFAULT NULL,
`previous_names`   blob  DEFAULT NULL,
`synonyms`         blob  DEFAULT NULL,
`accession_numbers`   blob  DEFAULT NULL,
`entrez_gene_id`   varchar(10)   DEFAULT NULL,
`ensembl_gene_id_hgnc`  varchar(40)   DEFAULT NULL,
`mouse_genome_database_id` varchar(40)   DEFAULT NULL,
`pubmed_ids`      blob  DEFAULT NULL,
`refseq_ids`      blob  DEFAULT NULL,
`uniprot_ids`     blob  DEFAULT NULL,
`ensembl_gene_id`  varchar(40)   DEFAULT NULL,
  PRIMARY KEY (`approved_symbol`),
  KEY `approved_symbol_idx` (`approved_symbol`),
  KEY `entrez_idx` (`entrez_gene_id`)
) ENGINE=MyISAM;
