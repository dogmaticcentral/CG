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
#include "utils.h"

/**********************************************/
int find_bin_index(Node *node, int value) {
    if (node->type==LEAF) {
        return node->id;
    } else if (value<=node->boundary) {
        return find_bin_index(node->left, value);
    } else {
        return find_bin_index(node->right, value);
    }

}

/**********************************************/
int set_boundaries(Node *node) {
    if (node->type == LEAF) {
        return node->boundary;
    } else {
        node->boundary = set_boundaries(node->left);
        return  set_boundaries(node->right);
    }
}

Node* inner_nodes(Node* node_list, int tier_start, int tier_end){
    if (tier_end-tier_start<=1) {
        node_list[tier_start].type = ROOT;
        return &node_list[tier_start];
    }
    /* n runs over indices in the node lis */
    /* t runs over indices in the tier */
    int n, t;
    n = tier_end; /* tier_end is upper bound on an open interval */
    for (t=tier_start; t<tier_end; t+=2) {
        node_list[n].type = INNER;
        node_list[n].id = n;
        node_list[n].left  = &node_list[t];
        node_list[n].left->parent = &node_list[n];
        if (t+1<tier_end) {
            node_list[n].right = &node_list[t+1];
            node_list[n].right->parent = &node_list[n];
        } else {
            node_list[n].right = node_list[n].left;
        }
        n++;
    }
    return inner_nodes(node_list, tier_end, n);
}

int tree_build_bottom_up (int boundaries [], int bdrs_length, Node* node_list, Node** root){
    int n;
    for (n=0; n< bdrs_length; n++) {
        node_list[n].type = LEAF;
        node_list[n].boundary = boundaries[n];
        node_list[n].id = n;
    }
    int tier_start=0; int tier_end=bdrs_length;
    *root = inner_nodes(node_list, tier_start, tier_end);
    set_boundaries(*root);
    return 0;
}

int number_of_nodes_needed(int bdrs_length) {

    int n = bdrs_length;
    int b = bdrs_length;
    while (b>1)  {
        b = (b+1)/2;
        n += b;
    }
    return n;
}

/**********************************************/
int  tree_init (Node *root, Node *left, Node *right) {

        printf ("Initalizing tree ...\n");

        root->type = ROOT;
        root->left      = left;
        root->right     = right;
        root->parent    = NULL;

        left->type = LEAF;
        left->left      = NULL;
        left->right     = NULL;
        left->parent    = root;

        right->type = LEAF;
        right->left      = NULL;
        right->right     = NULL;
        right->parent    = root;

        printf ("                 ... OK\n");
        return 0;
}

void node_print(FILE* fptr, Node* node ){
    fprintf ( fptr, "id:%3d     type:%3d   bdry:%3d     value:%10x    parent:%10x    left:%10x   right:%10x \n",
                  node->id,  node->type, node->boundary, (unsigned int)node, (unsigned int)node->parent,
                  (unsigned int)node->left,  (unsigned int)node->right);
}

int tree_print ( Node * node) {
        /* preorder print */
        if ( !node) return 1;

        node_print (stdout,  node);

        tree_print ( node->left);
        tree_print ( node->right);

        return 0;
}

int  tree_print_nhx (FILE* fptr,  Node * node){
        /*postorder*/
        if ( node->type == LEAF ) {
                fprintf ( fptr,  "%d", node->boundary );
        } else {
                fprintf ( fptr, "(" );
                if ( !node->left ) {
                        node_print (stderr, node);
                        PANIC ("Error in tree_print_nhx.\n");
                }
                tree_print_nhx ( fptr, node->left );
                fprintf ( fptr, "," );
                if ( !node->right ) {
                        node_print (stderr, node );
                        PANIC ("Error in tree_print_nhx.\n");
                }
                tree_print_nhx ( fptr, node->right );
                fprintf (fptr,  ")" );
        }
        if ( node->type == ROOT ) {
                fprintf (fptr, "\n");
        }
        return 0;

}

int tree_sanity_check ( Node * node) {
        /* preorder  */

        int retval=0, ok;

        if ( !node) return 0;

        switch (node->type) {
        case ROOT:
                ok = !(node->parent) &&  (node->left) && (node->right);
                if (!ok) retval = ROOT;
                break;
        case INNER:
                ok = (node->parent) &&  (node->left) && (node->right);
                if (!ok) retval = INNER;
                break;
        case LEAF:
                ok = (node->parent) &&  !(node->left) && !(node->right);
                if (!ok) retval = LEAF;
        }
        if (!ok ) {
                printf ("Sanity check failure:\n");
                node_print( stderr, node );
                return 1;
        }


        retval = tree_sanity_check ( node->left);
        if ( retval ) {
                return retval;
        }

        retval = tree_sanity_check ( node->right);
        if ( retval ) {
                return retval;
        }

        return 0;
}


int test ( ) {

        Node * root = NULL;
        Node * node_list = NULL;

        int boundaries [] = {1,2,5,7,10,12,13};
        int bds_length = 7;
        int nodes_needed = number_of_nodes_needed(bds_length);
        // int boundaries [] = {112,256};
        // int bds_length = 2;
        // int nodes_needed = bds_length + 1;

        node_list = (Node*)ecalloc(nodes_needed, sizeof(Node));
        tree_build_bottom_up(boundaries, bds_length, node_list, &root);

        Node * curr_root_ptr = root;
        int retval;
        if ( (retval=tree_sanity_check(curr_root_ptr)) ) {
                printf ("Tree trouble: %d.\n", retval);
                return 1;
        } else {
                printf ("Tree passes sanity check.\n");
                tree_print (curr_root_ptr );
                tree_print_nhx (stdout, curr_root_ptr );
        }
        printf("\n");
        int testvals [] = {0, 1, 3,  7, 8, 11, 13, 15};
        int t;
        for(t=0; t<8; t++) {
              printf("index for value %d:  %d\n", testvals[t], find_bin_index(root, testvals[t]));
        }
        return 0;
}
//int main ( int argc, char * argv[]) { test();}
