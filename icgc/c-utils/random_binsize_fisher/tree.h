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
# ifndef TREE_H
# define TREE_H


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define ROOT  1
# define INNER 2
# define LEAF  4

typedef struct Node {
        struct Node *left, *right, *parent;
        int id;
        int type;
        int boundary;
} Node;

int  number_of_nodes_needed(int number_of_leaves);
int  find_bin_index(Node *node, int value);
int  tree_build_bottom_up (int boundaries [], int bdrs_length, Node* node_list, Node** root);
int  tree_init (Node *root, Node *left, Node *right);
int  tree_sanity_check ( Node * node);
int  tree_print ( Node * node);
int  tree_print_nhx ( FILE * fptr, Node * node);

# endif
