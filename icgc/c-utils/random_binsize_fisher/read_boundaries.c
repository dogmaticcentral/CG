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

# include <stdio.h>
# include <stdlib.h>

#include "utils.h"
/****************************************************************************/
int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
             char * fmt, char * warnstr) {

        fprintf ( fptr, "Error on line %3d:     %s", line_ctr, line);
        fprintf ( fptr, fmt, warnstr);
        return 0;
}
/****************************************************************************/
int read_boundaries (char * filename, int ** bds_ptr, int * bds_size ) {
        FILE * fptr, *log = stdout;
        char line[LONGSTRING];
        char token[MAX_TOK][MEDSTRING] = {{'\0'}};
        char comment_char;
        int line_ctr = 0, retval;
        int max_token;
        /***************/

        fptr   = efopen ( filename, "r" );
        if (!fptr ) return 1;
        memset ( line, 0, LONGSTRING);
        *bds_size = 0;
        /* get the size and check the contents */
        while(fgets(line,LONGSTRING,fptr)!=NULL) {
                line_ctr++;
                /* tokenize */
                retval = tokenize ( token, &max_token, line, comment_char= '!' );
                switch ( retval ) {
                case  TOK_TOOMNY:
                        errmsg ( log, line_ctr, line, "\t\t %s\n", "Too many tokens.");
                        fclose (log);
                        break;
                case TOK_TOOLONG:
                        errmsg ( log, line_ctr, line, "\t\t %s\n", "Token too long.");
                        fclose (log);
                        break;
                }
                if ( max_token < 0 ) continue;
                (*bds_size) ++;
        }
        int* boundaries = (int*)ecalloc(*bds_size, sizeof(int));
        *bds_ptr = boundaries;
        rewind ( fptr);
        while(fgets(line,LONGSTRING,fptr)!=NULL) {
                tokenize ( token, &max_token, line, comment_char= '!' );
                if ( max_token < 0 ) continue;
                *boundaries = atoi(token[0]);
                boundaries++;
        }
        fclose(fptr);

        return 0;
}
