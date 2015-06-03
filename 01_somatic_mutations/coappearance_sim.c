# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

int **intmatrix(int rows, int columns){
    int **m;
    int i;
        /* allocate pointers to rows */
    m=(int **) malloc(rows*sizeof(int*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(int *) calloc( rows*columns, sizeof(int));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

int main (int argc, char * argv[]) {

    if ( argc  < 6) {
	fprintf (stderr,
		 "Usage: %s <no_slots> <no_red> <no_black> <observed_no_coappearances> <no_iterations>\n",
		 argv[0]);
	exit(1);
    }
    int M  = atoi(argv[1]);
    int Nr = atoi(argv[2]);
    int Nb = atoi(argv[3]);
    int coapps               = atoi(argv[4]);
    int number_of_iterations = atoi(argv[5]);

    srand48(time(NULL));

    //printf ("M = %d, Nr = %d, Nb = %d, coapps = %d, number_of_iterations = %d\n",
    //M, Nr, Nb, coapps, number_of_iterations);

    double avg_number_of_double_labeled = 0;
    double pval_le = 0.0;   // probabilty of being less-or_equal-to observed
    double pval_ge = 0.0;   // probabilty of being greater-or-equal-to observed

    if (number_of_iterations <= 0  || Nr<=0 || Nb <= 0) {
	printf ("  %10.2lf   %10.2lf  %10.2lf  \n",  avg_number_of_double_labeled, pval_le, pval_ge);
	return 0;
    }
    // slot[x][0] correponds to the number of red marbles in slot x
    // slot[x][1] correponds to the number of black marbles in slot x
    int ** slots = intmatrix(M,2);
    int i;
    for (i=0; i < number_of_iterations; i++) {
	
	//#####
        memset ( slots[0], 0, M*2*sizeof(int) );
        
        int number_of_double_labeled = 0;
	int r, b, s;
	for (r=0; r<Nr; r++) {
            int random_slot = (int)(drand48()*M);
            slots[random_slot][0] += 1;
	}
        for (b=0; b<Nb; b++) {
            int random_slot = (int)(drand48()*M);
            slots[random_slot][1] += 1;
	}
	for (s=0; s<M; s++) {
            if ( slots[s][0] &&  slots[s][1])   number_of_double_labeled += 1;
	}

        //#####
        avg_number_of_double_labeled += number_of_double_labeled;
        if ( number_of_double_labeled <= coapps)   pval_le += 1.0;
        if ( number_of_double_labeled >= coapps )  pval_ge += 1.0;
    }
    //##################################
    avg_number_of_double_labeled /= (double)(number_of_iterations);
    pval_le /= (double)(number_of_iterations);
    pval_ge /= (double)(number_of_iterations);

    printf ("  %10.2lf     %10.5lf   %10.5lf  \n",  avg_number_of_double_labeled, pval_le, pval_ge);
    
    return 0;
}
