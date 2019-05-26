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

# ifndef UMBRELLA_H
# define UMBRELLA_H

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "pdb.h"
# include "utils.h"

# define PANIC(msg, str) {			\
    fprintf (stderr, "%s %s.\n",msg, str);	\
    exit(1);					\
}

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
} Atom;
# define  MAX_NO_ATOMS 100
typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+1];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
} Residue;

Residue * sequence1, *sequence2;
int no_res_1, no_res_2;

typedef struct {
    char aa;
    double distance;
} Point;


extern char *amino_acid_order;
# define AA_ORDER_LENGTH 21

int read_pdb ( char * pdbname, char *chain_id_ptr, Residue ** sequence_ptr, int * no_res_ptr);

int read_selection (Residue *sequence, int no_res, char * filename, int * selection);
int determine_dist_matrix ( double ** dist_matrix, Residue * sequence, int no_res);

void cluster_counter (int  no_of_things,  int *neighbors[], int * mask,
		     int cluster_count_per_size[], int * no_of_clusters,
		      int * max_size, int * secnd_max_size , int * clusters[]);
int clustering_z_score ( int no_res, int ** neighbors, int * selected,  double *clustering_score);

# endif
