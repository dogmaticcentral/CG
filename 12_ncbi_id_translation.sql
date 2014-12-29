/* the cooresponding file  in ~/database directory is called gene_info_human_ncbi.txt*
/* tha table is supposed to go to baseline database */
CREATE TABLE ncbi_id_translation (
  tax_id int NOT NULL, 
  gene_id INT NOT NULL, 
  symbol  varchar(40)   DEFAULT NULL,
  locus_tag  blob  DEFAULT NULL,
  synonyms blob  DEFAULT NULL,
  db_xrefs blob  DEFAULT NULL,
  chromosome varchar (10), 
  map_location blob  DEFAULT NULL,
  description blob  DEFAULT NULL,
  type_of_gene blob  DEFAULT NULL,
  symbol_from_nomenclature_authority blob  DEFAULT NULL,
  full_name_from_nomenclature_authority blob  DEFAULT NULL,
  nomenclature_status blob  DEFAULT NULL,
  other_designations blob  DEFAULT NULL,
  modification_date   blob  DEFAULT NULL,
   PRIMARY KEY (gene_id)
) ENGINE=MyISAM;

create index symbol_idx on ncbi_id_translation (symbol);
