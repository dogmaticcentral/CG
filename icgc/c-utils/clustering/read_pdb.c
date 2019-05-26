/*
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
#
# Original publication:  https://www.ncbi.nlm.nih.gov/pubmed/12875851
#
#
*/


# include <ctype.h>
# include "pdbclust.h"
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

# define BUFFLEN 150
int read_pdb ( char * pdbname, char *chain_id_ptr, Residue ** sequence_ptr, int * no_res_ptr) {

  Residue * sequence;
  FILE * fptr = NULL;
  char line[BUFFLEN];
  char oldresno[PDB_ATOM_RES_NO_LEN+1];
  char oldrestype [PDB_ATOM_RES_NAME_LEN+1];
  char tmp[PDB_ATOM_X_LEN+1];
  int atomctr, resctr, no_res, skip;
  char chain_id = *chain_id_ptr;

  /* open file */
  fptr = fopen ( pdbname, "r");
  if ( !fptr ) {
    fprintf (stdout, "Cno %s.\n", pdbname);
    return 1;
  }
  /* check for NMR and theory files */
  while(fgets(line, BUFFLEN, fptr)!=NULL){
      if( ! strncmp(line,"EXPDTA", 6)){
          if ( strstr (line, "NMR") ) {
            fprintf (stdout, "%s is an NMR file.\n", pdbname);
            return 1;
          } else if (  strstr (line, "THEO") ) {
            fprintf (stdout, "%s is a theoretical  file.\n", pdbname);
            return 1;
          }
      }
  }

  /* count residues */
  memset (line,  0, BUFFLEN);
  memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+1);
  resctr = 0;
  rewind ( fptr);
  while(fgets(line, BUFFLEN, fptr)!=NULL){
    if( ! strncmp(line,"ATOM", 4)   ||    ! strncmp(line,"HETATM", 6)     ){
      /* if chain specified, check if it is the right one: */
      if (chain_id  &&   chain_id != line[PDB_ATOM_CHAINID] )  continue;

      if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN) ) {

        strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN);
        oldresno[PDB_ATOM_RES_NO_LEN] = '\0';
        /* printf ( "New residue number:  %s \n", oldresno); */
        resctr ++;
      }

    }
  }

  no_res = resctr;
  *no_res_ptr = no_res;
  printf ("reading chain %c   no residues: %d\n", chain_id, no_res);

  /* allocate space */
  sequence = NULL;
  sequence = calloc ( no_res, sizeof (Residue));
  if ( ! sequence ) {
    fprintf ( stderr, "Error allocating sequence space.\n");
    exit (1);
  }
  *sequence_ptr = sequence;

  /* read in the atoms */
  rewind ( fptr);
  memset (line,  0, BUFFLEN);
  memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+1);
  memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+1);
  resctr= -1;
  atomctr = 0;
  skip = 0;
  while(fgets(line, BUFFLEN, fptr)!=NULL){
        if ( skip ) { /*skipping the alt locations; see below */
          skip = 0;
          continue;
        }
        if( ! strncmp(line,"ATOM", 4)   ||    ! strncmp(line,"HETATM", 6)     ){

            /* if chain specified, check if it is the right one: */
            if (chain_id  &&   chain_id != line[PDB_ATOM_CHAINID] )  continue;

            /* skip hydrohens */
            if ( line[PDB_ATOM_ATOM_NAME] == 'H' ||  line[PDB_ATOM_ATOM_NAME+1] == 'H') {
              continue;
            }
            /* adjust the counters */
            if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN) ) {
                  strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN);
                  strncpy (oldrestype, line+PDB_ATOM_RES_NAME, PDB_ATOM_RES_NAME_LEN);
                  oldresno[PDB_ATOM_RES_NO_LEN] = '\0';
                  oldrestype[PDB_ATOM_RES_NAME_LEN] = '\0';
                  resctr ++;
                  atomctr = 0;

                  sequence[resctr].no_atoms = 1;
                  strncpy ( sequence[resctr].pdb_id,oldresno, PDB_ATOM_RES_NO_LEN);
                  sequence[resctr].pdb_id[PDB_ATOM_RES_NO_LEN] = '\0';
                  strncpy ( sequence[resctr].res_type,oldrestype, PDB_ATOM_RES_NAME_LEN);
                  sequence[resctr].res_type[PDB_ATOM_RES_NAME_LEN] = '\0';

            } else {
                  atomctr ++;
                  sequence[resctr].no_atoms = atomctr + 1;
                  if ( atomctr >= MAX_NO_ATOMS ) {
                    fprintf ( stdout, "Error: I thought every aa has < %d atoms.\n",
                    MAX_NO_ATOMS ); /*stdout so perl can capture it easily */
                    exit (1);
                  }
            }

            /* read in atom info */
            strncpy ( sequence[resctr].atom[atomctr].type, line+ PDB_ATOM_ATOM_NAME, PDB_ATOM_ATOM_NAME_LEN );
            sequence[resctr].atom[atomctr].type[ PDB_ATOM_ATOM_NAME_LEN] = '\0';
            strncpy ( tmp, line+PDB_ATOM_X , PDB_ATOM_X_LEN);
            tmp[PDB_ATOM_X_LEN] = '\0';
            sequence[resctr].atom[atomctr].x=atof(tmp);
            strncpy ( tmp, line+PDB_ATOM_Y , PDB_ATOM_Y_LEN );
            tmp[PDB_ATOM_Y_LEN] = '\0';
            sequence[resctr].atom[atomctr].y=atof(tmp);
            strncpy ( tmp, line+PDB_ATOM_Z , PDB_ATOM_Z_LEN);
            tmp[PDB_ATOM_Z_LEN] = '\0';
            sequence[resctr].atom[atomctr].z=atof(tmp);
            if ( ! isspace (line[PDB_ATOM_ALTLOC]) ) {
                /* atom   has alternative location  */
                /* skip the next line */
                skip = 1;
            }
        } // end of ATOM, HETATM line processing
  } // end of while readline

  /* clean PDB id tags from spaces */
  for (resctr=0; resctr < no_res; resctr ++ ) {
        if ( string_clean (sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN+1)) {
          fprintf (stderr, "Error in read_pdb(): empty id string for residue with sequential no %d.\n", resctr);
          return 1;
        }
        if ( string_clean (sequence[resctr].res_type, PDB_ATOM_RES_NAME_LEN)) {
          fprintf (stderr, "Error in read_pdb(): empty id string for residue with sequential no %d.\n", resctr);
          return 1;
        }
        sequence[resctr].res_type_short  = single_letter ( sequence[resctr].res_type );
  }

  /* close file */
  fclose (fptr);
  return 0;
}
