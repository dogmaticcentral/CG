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

#include "pdbclust.h"

int clustering_z_score ( int no_res, int ** neighbors, int * selected,  double *clustering_score){

    int pos;
    int first, no_selected;
    double score, avg, std_dev, z;
    int cluster_score (int no_of_res, int *seq, int ** adj_matrix,double *score);
    int std_dev_over_S (int L, int M, int ** adj_matrix, double *avg, double * std_dev, int first);

    no_selected = 0;
    for (pos=0; pos< no_res; pos++) no_selected += selected[pos];
    
    /* find  scpre */
    cluster_score (no_res, selected, neighbors, &score);
    /* find avg and stddev in the set of random picks */
    std_dev_over_S (no_res, no_selected, neighbors, &avg,  &std_dev, first = 1);
	
    /* evaluate and store the z-score */
    z = (std_dev>1.e-5) ? (score - avg)/std_dev : 0.0;
    *clustering_score  = z;
    
    return 0;
}

/************************************************************************************/
int cluster_score (int no_of_res, int *selected, int ** neighbors, double *score){
    int i,j, dist;
    int  size;
    double sum;

    size = no_of_res;
    sum = 0.0;
    for (i=0; i<size-1 ; i++) {	
	if ( !selected[i]) continue;
	for (j=i+1; j<size; j++) {
	    if (selected[j] && neighbors[i][j]) {
		dist = j-i; // note: we are ingoring the fact there might be gaps in the structure
		sum += dist;
	    }
	}
	
    }
    *score = sum;
    return 0;    
}

/************************************************************************************/
int std_dev_over_S (int no_res, int no_selected, int ** neighbors, double *avg, double * std_dev, int first) {
    
    int i,j,k,l,n;
    double  std_dev_thry, avg_thry;
    double ratio[3];
    double aux;
    static double subsum[3], bare_avg_thry;

    if (first) {
	bare_avg_thry = 0;
	subsum[0] = subsum[1] = subsum[2] = 0.0;
	for (i=0; i<no_res-1; i++) {
	    for (j=i+1; j<no_res; j++) {
			
		if ( neighbors[i][j]) {
			    
		    bare_avg_thry +=  (j-i);
			    
		    for (k=0; k<no_res-1; k++) {
			for (l=k+1; l<no_res; l++) {
				    
			    if  ( neighbors[k][l]) {
				n = (k==i) + (l==j) + (k==j) + (l==i);
				subsum [2-n] += (j-i)*(l-k);
			    }
			}
		    }
		}
	    }
	}
        bare_avg_thry /= (no_res*(no_res-1));
    }
    
    ratio[0] = (double) no_selected*(no_selected-1)/(no_res*(no_res-1));
    ratio[1] = ratio[0]*(no_selected-2)/(no_res-2);
    ratio[2] = ratio[1]*(no_selected-3)/(no_res-3);
    std_dev_thry = 0.0;
    for (i =0; i<3; i++) {
	std_dev_thry += subsum[i]*ratio[i];
    }
   
    avg_thry = bare_avg_thry*no_selected*(no_selected-1);

    aux = (std_dev_thry-avg_thry*avg_thry)/std_dev_thry;
    if ( aux  < -1.0e-4){
	fprintf (stderr,"Unspecified error in std dev (S) calculation.\n");
	fprintf (stderr, " %8.3e  %8.3e  %8.3e \n", std_dev_thry, avg_thry*avg_thry, aux);
	return 1;
    } else if ( aux < 0 ) {
	*std_dev = 0;
    } else {
	std_dev_thry -= avg_thry*avg_thry;
	*std_dev = sqrt(std_dev_thry);
    }
    *avg = avg_thry;

    return 0;
    
}

