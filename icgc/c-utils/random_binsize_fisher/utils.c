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

#include "utils.h"

void error (int errno, char *errstr) {
        fprintf (stderr, "%s\n", errstr);
        exit (errno);
}

void usage(char *  use[])
{

        if (use != NULL) {
                (void) fprintf(stdout, "\n\t%s\n\n", *use);
                while (*++use != NULL)
                        (void) fprintf(stdout, "\t%s\n", *use);
                (void) fprintf(stdout, "\n");
        }

        return;
}


void * emalloc(size_t size) {
        void * ptr;
        if ((ptr = calloc(size, 1)) == NULL) {
                fprintf (stderr,  "emalloc: no memory for %zu bytes", size);
                exit (1);
        }
        return ptr;
}

/**********************************************/
void *  ecalloc (size_t nitems, size_t size) {
        void * retval;
        if (  (retval=  calloc(nitems, size))  ) {
                return retval;
        } else {
                PANIC ("Memory allocation failure.\n");
                return NULL;
        }
}


FILE * efopen(char * name, char * mode)
{

        FILE * fp;


        if ((fp = fopen(name, mode)) == NULL) {
                fprintf (stderr,
                         "efopen: can't open \"%s\" for \"%s\"\n", name, mode);
                exit (1);
        }

        return fp;

}






/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch]
 */
char **chmatrix(long nrl, long nrh, long ncl, long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **m;

        /* allocate pointers to rows */
        m=(char **) malloc((size_t)((nrow+1)*sizeof(char*)));
        if (!m)  {
                fprintf (stderr,"allocation failure 1 in matrix()");
        }
        m += 1;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(char *) calloc( nrow*ncol+1, sizeof(char));
        if (!m[nrl]) {
                fprintf (stderr,"allocation failure 2 in matrix()for char");
        }
        m[nrl] += 1;
        m[nrl] -= ncl;

        for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}



/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]
 */
int **imatrix(long nrl, long nrh, long ncl, long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int **m;

        /* allocate pointers to rows */
        m=(int **) malloc((size_t)((nrow+1)*sizeof(int*)));
        if (!m) fprintf (stderr,"allocation failure 1 in matrix()");
        m += 1;
        m -= nrl;


        /* allocate rows and set pointers to them */
        m[nrl]=(int *) calloc(nrow*ncol+1,sizeof(int));
        if (!m[nrl]) fprintf (stderr,"allocation failure 2 in matrix() for int");
        m[nrl] += 1;
        m[nrl] -= ncl;

        for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}


/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
 */
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
        if (!m) fprintf (stderr,"allocation failure 1 in matrix()");
        m += 1;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) calloc( nrow*ncol+1,sizeof(double));
        if (!m[nrl]) fprintf (stderr,"allocation failure 2 in matrix()for double");
        m[nrl] += 1;
        m[nrl] -= ncl;

        for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}



/* allocate a float square matrix
 */
float **fmatrix(long dimension)
{
        long i;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((size_t)((dimension)*sizeof(float*)));
        if (!m) fprintf (stderr,"allocation failure 1 in matrix()");

        for(i=0; i<dimension; i++)
        {
                m[i] = (float *) calloc(dimension,sizeof(float));
                if (!m[i]) fprintf (stderr,"allocation failure 1 in matrix()");
        }

        /* return pointer to array of pointers to rows */
        return m;
}



/* free a char matrix allocated by chmatrix()
 */
void free_chmatrix(char **m, long nrl, long nrh, long ncl, long nch)
{
        free( (m[nrl]+ncl-1));
        free( (m+nrl-1));
}



/* free an int matrix allocated by imatrix()
 */
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
        free( (m[nrl]+ncl-1));
        free( (m+nrl-1));
}



/* free a double matrix allocated by dmatrix()
 */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
        free( (m[nrl]+ncl-1));
        free( (m+nrl-1));
}




/* sort array according to the score in the other */
/* I couldn't declare pos_cmp within array_qsort  bcs it   crashed on mac */

double * score_array;

int pos_cmp (const void * a0, const void * b0) {

        int * a= (int*) a0;
        int * b= (int*)b0;
        if ( score_array[*a] > score_array[*b]) {
                return 1;
        }
        if ( score_array[*a] < score_array[*b]) {
                return -1;
        }
        return 0;
}


int array_qsort (int * sorted_pos, double * sa, int sequence_length ) {
        /* position comparison function */
        score_array = sa;

        qsort (sorted_pos+1, sequence_length, sizeof(int), pos_cmp);

        return 0;
}


/***************************************************************************/
int tokenize ( char token[MAX_TOK][MEDSTRING], int * max_token,
               char * line, char comment_char) {
        /* assumes the tokens to be no bigger than MEDSTRING */

        char * chrptr, *last;
        int current_token, current_char = 0;
        int reading;

        memset (token[0], 0, MAX_TOK*MEDSTRING*sizeof(char));
        chrptr = line;
        last   = chrptr + strlen (line);
        current_token = -1;
        current_char  =  0;
        reading = 0;
        while ( chrptr <= last) {
                if ( *chrptr == comment_char ) break;
                if ( *chrptr == '\n' ) break;
                if ( *chrptr && !isspace(*chrptr) ) {
                        if ( !reading ) {
                                reading = 1;
                                current_char = 0;
                                current_token++;
                                if ( current_token >= MAX_TOK ) {
                                        return TOK_TOOMNY; /* defined in possum_utils.h */
                                }
                        }
                        if ( current_char >= MEDSTRING ) {
                                return TOK_TOOLONG;
                        }
                        token[current_token][current_char] = *chrptr;
                        current_char++;
                } else {
                        if ( reading ) {
                                reading = 0;
                        }
                }
                chrptr++;
        }
        *max_token = current_token;

        return 0;

}

/**********************************************************/
/* get rid of spaces in a string */
int  string_clean ( char* string, int length) {
        int ctr;
        for (ctr = 0; ctr < length; ctr++) {
                if ( isspace (string[ctr]) ) string[ctr] = '\0';
        }
        ctr=0;
        while ( !string[ctr] && ctr < length) ctr++;

        if ( ctr == length ) return 1; /* empty string */

        if ( ctr ) {
                memmove (string, string+ctr, length-ctr);
                memset ( string+length-1-ctr, 0, ctr);
        }

        return 0;
}

/**********************************/
char single_letter ( char code[]){

        switch ( code[0] ) {
        case 'A':
                switch ( code [1]) {
                case 'L':
                        return 'A';
                        break;
                case 'R':
                        return 'R';
                        break;
                case 'S':
                        switch ( code[2] ) {
                        case 'N':
                                return 'N';
                                break;
                        case 'P':
                                return 'D';
                                break;
                        }
                        break;
                }
                break;
        case 'C':
                return 'C';
                break;
        case 'G':
                /* the second letter is always L */
                switch ( code[2] ) {
                case 'U':
                        return 'E';
                        break;
                case 'N':
                        return 'Q';
                        break;
                case 'Y':
                        return 'G';
                        break;
                }
                break;
        case 'H':
                return 'H';
                break;
        case 'I':
                return 'I';
                break;
        case 'L':
                switch ( code [1]) {
                case 'E':
                        return 'L';
                        break;
                case 'Y':
                        return 'K';
                        break;
                }
                break;
        case 'M':
                return 'M';
                break;
        case 'P':
                switch ( code [1]) {
                case 'H':
                        return 'F';
                        break;
                case 'R':
                        return 'P';
                        break;
                case 'T':
                        return 'Y'; //PTR phosphotyrosine
                        break;
                }
                break;
        case 'S':
                switch ( code [1]) {
                case 'E':
                        return 'S';
                        break;
                case 'C':
                        return 'C'; // SCY, selenocysteine
                        break;
                }
        case 'T':
                switch ( code [1]) {
                case 'H':
                        return 'T';
                        break;
                case 'R':
                        return 'W';
                        break;
                case 'Y':
                        return 'Y';
                        break;
                case 'P':
                        return 'T'; // TPO phosphothreonine
                        break;
                }
                break;
        case 'V':
                return 'V';
                break;

        }


        fprintf (stdout, "Unrecognized amino acid code: %s.\n", code);
        return 0;
}
