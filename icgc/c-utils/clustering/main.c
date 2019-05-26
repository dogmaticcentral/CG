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

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include "utils.h"
# include "pdbclust.h"


/************************************************************/
int main ( int argc, char * argv[]) {

    char pdbname[150] = {'\0'};
    char selection_file[150] = {'\0'};
    Residue * sequence;
    int no_res, res_ctr;
    int * selected;
    double cutoff_dist, score;
    double **distmat;
    int ** cluster;
    int c;
    char chain_id = '\0';
    FILE * fclust = NULL;


    if ( argc < 5 ) {
      	fprintf (stderr,
      		 "Usage: %s <pdbfile>  <chain id>|'-'  <selection file>  <cutoff dist> \n"
           "where selection file is of the format (for example) \n"
           "F    3 \n"
           "K   48 \n"
           "R   50 \n"
           "D   59 \n"
           "E   82 \n"
           "G   91 \n"
           "A   97 \n"
           "C  100 \n"
           "etc \n"
           "and the cutoff distance is given in Angstroms\n",
      		 argv[0]);
      	exit (1);
    }
    sprintf ( pdbname, "%s", argv[1]);
    chain_id =  argv[2][0]=='-' ? '\0' : argv[2][0];
    sprintf ( selection_file, "%s", argv[3]);
    cutoff_dist = atof ( argv[4]);
    printf ("cutoff distance: %5.2lf Angstroms\n", cutoff_dist);

    /* input the structure */
    if ( read_pdb ( pdbname, &chain_id, &sequence, &no_res) ) exit (1);
    printf ("there are %d residues in %s, chain %c.\n",  no_res, pdbname, chain_id);

    /* input selection */
    selected = emalloc (no_res*sizeof(int));
    if (! selected) {
      	fprintf (stderr, "error allocating selection array\n");
      	exit(1); // leaky, leaky
    }
    if ( read_selection (sequence, no_res, selection_file, selected) ) exit (1);

    /* allocate space for dist matrix */
    distmat = (double**) dmatrix (0,  no_res-1, 0,  no_res-1);
    /* calculate dist */
    if ( determine_dist_matrix(distmat,  sequence, no_res) ) exit(1);

    /* cluster counting ... */
  	int  no_of_clusters, max_size, secnd_max_size;
  	int * cluster_count;
  	int ** neighbors;
  	int res1, res2;

  	cluster_count       =  (int *) emalloc ( (no_res+1)*sizeof(int));
  	cluster             =  imatrix (0, no_res, 0, no_res);
  	neighbors           =  imatrix (0, no_res-1, 0, no_res-1);

  	for (res1=0;  res1 < no_res; res1++ ) {
  	    neighbors [res1][res1] = 1;
  	    for (res2= res1+1; res2 < no_res; res2++ ) {
        		neighbors[res1][res2] = ( distmat[res1][res2] < cutoff_dist);
        		neighbors[res2][res1] = neighbors[res1][res2];
  	    }
  	}

	  cluster_counter (no_res,  neighbors,  selected, cluster_count, & no_of_clusters,
			   &max_size, &secnd_max_size , cluster);
	  clustering_z_score ( no_res,  neighbors,  selected, &score);

  	if (0) {
        	  printf ("runnning simulation ...\n");
      	    int cluster_score (int no_of_res, int *seq, int ** adj_matrix,double *score);
      	    int i, n;
      	    int no_selected = 0;
      	    for (i=0; i<no_res; i++) no_selected += selected[i];
      	    cluster_score (no_res, selected, neighbors, &score);
      	    printf (" selected  %4d   original cluster score = %8.2e\n", no_selected, score);
      	    double frac = (double)no_selected/no_res;
      	    double score = 0;
      	    for (n=0; n<100; n++) { // number of reps
      		      memset (selected, 0, no_res*sizeof(int));
            		for (i=0; i<no_res; i++) {
            		    if (drand48() < frac) selected[i] = 1;
            		}
            		cluster_score (no_res, selected, neighbors, &score);
            		no_selected = 0;
            		for (i=0; i<no_res; i++) no_selected += selected[i];
            		printf  (" selected  %4d    random score = %8.2e\n",  no_selected, score);
      	    }
  	}

    /* output */
    fclust = stdout;
    for ( c=0; c <= no_res; c++) {
      	if ( ! cluster[c][0] ) {
      	    continue;
      	}
      	if ( !c ) {
      	    fprintf ( fclust,"\t isolated:\n");
      	} else {
      	    fprintf ( fclust,"\t cluster size: %3d \n", cluster[c][0]);
      	}
      	for ( res_ctr=1; res_ctr <=  cluster[c][0]; res_ctr++) {
      	    fprintf ( fclust, "%s  \n", sequence[ cluster[c][res_ctr] ].pdb_id );
      	}
    }
    fprintf (fclust, "\nclustering  z-score:  %8.3f\n", score);



    return 0;



}
