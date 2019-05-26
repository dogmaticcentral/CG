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

# ifndef _UTILS_H
# define _UTILS_H
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>

#include <sys/types.h>
# include <stdio.h>
#include <ctype.h>

# define BUFFLEN     150
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25


/******************************/
/*   tokenizer                */
/******************************/

# define TOK_TOOMNY  1 /* tokenizer error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */

/****************************/
void usage(char *  use[]);
void error (int errno, char *errstr);
void * emalloc(size_t	size);
FILE * efopen(char * name, char * mode);
char **chmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
float **fmatrix(long dimension);
void free_chmatrix(char **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

int array_qsort (int * sorted_pos, double * sa, int sequence_length );

int tokenize ( char token[MAX_TOK][MEDSTRING], int * max_token,
	       char * line , char comment_char);
int  string_clean ( char* string, int length);
char single_letter ( char code[]);
# endif
