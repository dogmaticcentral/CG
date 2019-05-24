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
*/
# include "tree.h"
# include "utils.h"
# include <time.h>

int read_boundaries (char * filename, int ** bds_ptr, int * bds_size );

/*****************************************************/
int mark_selection(Node * root, int nsel, char * selection_array, int max_val, char flag) {
    int ctr, randval, rand_indx;
    for (ctr=0; ctr<nsel; ctr++) {
        do {
            // drand48 generates numbers on closed interval [0.0, 1.0]
            randval   = (int)round(max_val*drand48());
            rand_indx = find_bin_index(root, randval);
        } while (selection_array[rand_indx]&flag); // no repeats
        selection_array[rand_indx] |= flag;
    }
    return 0;
}
/*****************************************************/
int main ( int argc, char * argv[]) {
    /* parse command line */
    if ( argc < 6 ) {
      	fprintf (stderr,
      		 "Usage: %s  <bin bdries file>  <selection size 1>   <selection size 2>  <overlap size> <number of simluation rounds>\n"
           "where <bin bdries file is of the format (for example) \n"
           " 113 \n"
           " 256 \n"
           " 500 \n"
           "etc \n",
      		 argv[0]);
      	exit (1);
    }
    char infilename[150] = {'\0'};
    int nsel1, nsel2, noverlap, nrounds;
    sprintf ( infilename, "%s", argv[1]);
    nsel1 = atoi(argv[2]);
    nsel2 = atoi(argv[3]);
    noverlap = atoi(argv[4]);
    nrounds = atoi(argv[5]);
    /* read in boundary array */
    int * boundaries;
    int bds_length;
    if (read_boundaries (infilename, &boundaries,  &bds_length)) exit (1);
    if (nsel1>bds_length || nsel2>bds_length) {
        fprintf (stderr,  "Both selection sizes must be no bigger than the number of bins.\n");
        exit(1);
    }

    /* build the search tree */
    int nodes_needed = number_of_nodes_needed(bds_length);
    Node * root = NULL;
    Node * node_list = NULL;
    node_list = (Node*)ecalloc(nodes_needed, sizeof(Node));
    tree_build_bottom_up(boundaries, bds_length, node_list, &root);

    /* simulation */
    srand48(time(NULL));
    int n, i;
    char * selection_array = calloc(bds_length, sizeof(char));
    int max_val = boundaries[bds_length-1];
    int count_smaller=0, count_bigger=0;
    for(n=0; n<nrounds; n++) {
        memset(selection_array, 0, bds_length*sizeof(char));
        mark_selection(root, nsel1, selection_array, max_val, 1);
        mark_selection(root, nsel2, selection_array, max_val, 2);
        int overlap_size = 0;
        for (i=0; i<bds_length; i++) {
            if (selection_array[i]==3) overlap_size+=1;
        }
        if (overlap_size>=noverlap) count_bigger++;
        if (overlap_size<=noverlap) count_smaller++;
    }
    count_smaller +=1; // to remind ourselves that e cannot go below 1/nrounds in precision
    printf("OK\t%.2e\t%.2e\n", (float)count_smaller/nrounds, (float)count_bigger/nrounds);
    return 0;
}
