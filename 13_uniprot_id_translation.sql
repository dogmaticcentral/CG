/* from ~/databases/uniprot_human.idmapping.dat */
/* to baseline database*/
CREATE TABLE uniprot_id_translation (
  uniprot_id  varchar(10)   NOT NULL, 
  other_db    varchar(40)   DEFAULT NULL,
  other_db_id varchar(100) DEFAULT NULL
) ENGINE=MyISAM;

create index uniprot_idx on uniprot_id_translation (uniprot_id);
create index other_idx on uniprot_id_translation (other_db_id);
