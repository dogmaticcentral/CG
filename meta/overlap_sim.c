# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <time.h>


int main ( int argc, char * argv[]) {

    if ( argc < 6) {
	fprintf (stderr, "Usage: %s <no of slots> <no of red labels> <no of black labels> <overlap observed> <no of iterations>.\n",  argv[0]);
	exit(1);
    }
    
    int M, Nr, Nb, obs, number_of_iterations;
    M  = atoi(argv[1]);
    Nr = atoi(argv[2]);
    Nb = atoi(argv[3]);
    obs  = atoi(argv[4]);
    number_of_iterations = atoi(argv[5]);

    double avg_number_of_double_labeled = 0;
    double pval_le = 0.0;
    double pval_ge = 0.0; 
   

    if (number_of_iterations <=0 ) {
	printf ( " %8.2lf    %8.6lf       %8.6lf   \n", avg_number_of_double_labeled, pval_le, pval_ge);
	return 0;
    }
    srand48(time(NULL));
    int i, s, slots[M][2], number_of_double_labeled;
    int random_slot;
    for (i=0; i<= number_of_iterations; i++ ) {

	number_of_double_labeled = 0;
	memset (slots[0], 0, M*2*sizeof(int) );
	for (s=0; s<Nr; s++) {
	    random_slot = drand48()*M;
	    slots[random_slot][0] += 1;
	}
	for (s=0; s<Nb; s++) {
	    random_slot = drand48()*M;
	    slots[random_slot][1] += 1;
	}
	for (s=0; s<M; s++) {
	    if (slots[s][0] >0 && slots[s][1] > 0 ) number_of_double_labeled += 1;
	}
	
	avg_number_of_double_labeled += number_of_double_labeled;
	if ( number_of_double_labeled <= obs ) pval_le += 1.0;
        if ( number_of_double_labeled >= obs)  pval_ge += 1.0;

    }
    avg_number_of_double_labeled /= (float)number_of_iterations;
    pval_le /= (float)number_of_iterations;
    pval_ge /= (float)number_of_iterations;

    float eval_le, eval_ge;
    if (pval_le<=0) {
	eval_le= 100.00;
    } else if (pval_le>=1) {
	eval_le = 0.0;
    } else {
	eval_le = -log10(pval_le);
    }
    if (pval_ge<=0) {
	eval_ge= 100.00;
    } else if (pval_ge>=1) {
	eval_ge = 0.0;
    } else {
	eval_ge = -log10(pval_ge);
    }
    
    printf ( " %8.3lf   %8.2lf   %8.2lf   \n", avg_number_of_double_labeled, eval_le, eval_ge);

    
    return 0;
}
