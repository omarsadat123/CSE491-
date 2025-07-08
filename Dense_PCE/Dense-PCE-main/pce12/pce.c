/* PCE: Pseudo Clique Enumerater */
/* ver 1.0 1/Sep/2005 Takeaki Uno   e-mail:uno@nii.jp, 
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, do not forget to 
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about the code for the users.
   For the commercial use, please make a contact to Takeaki Uno. */


#ifndef _pce_c_
#define _pce_c_

#define WEIGHT_DOUBLE

#include"alist.c"
#include"sgraph.c"
#include"problem.c"

// int PCE_connected=1;
void PCE_error (){
  ERROR_MES = "command explanation";
  print_err ("pce: CMqsS [options] input-filename threshhold [output-filename]\n\
%%:show progress, _:no message, +:write solutions in append mode\n\
C:enumerate pseudo cliques, M:enumerate maximal pseudo cliques\n\
Q,f:print density preceding/following to pseudo cliques\n\
[options]\n\
-K [num]:output [num] most frequent itemsets\n\
-l [num]:output pseudo cliques with size at least [num]\n\
-u [num]:output pseudo cliques with size at most [num]\n\
-U [real-num]:upper bound for density\n\
-# [num]:stop after outputting [num] solutions\n\
-, [char]:give the separator of the numbers in the output\n\
-Q [filename]:replace the output numbers according to the permutation table given by [filename]\n\
threshold is given by real-number from 0 to 1\n\
# the 1st letter of input-filename cannot be '-'.\n\
# if the output file name is '-', the solutions will be output to standard output.\n");
  EXIT;
}
// c:enumerate non-connected pseudo cliques\n");
// -w [filename]:read weights of edges from [filename]\n");

/***********************************************************************/
/*  read parameters given by command line  */
/***********************************************************************/
void PCE_read_param (int argc, char *argv[], PROBLEM *PP){
  ITEMSET *II = &PP->II;
  int c=1;
  if ( argc < c+3 ){ PCE_error (); return; }

  if ( !strchr (argv[c], '_') ){ II->flag |= SHOW_MESSAGE; PP->SG.flag |= SHOW_MESSAGE; }
  if ( strchr (argv[c], '%') ) II->flag |= SHOW_PROGRESS;
  if ( strchr (argv[c], '+') ) II->flag |= ITEMSET_APPEND;
  if ( strchr (argv[c], 'f') ) II->flag |= ITEMSET_FREQ;
  if ( strchr (argv[c], 'Q') ) II->flag |= ITEMSET_PRE_FREQ;
  if ( strchr (argv[c], 'C') ) PP->problem |= PROBLEM_FREQSET;
//  if ( strchr( argv[c], 'c' ) ) PP->problem |= FREQSET;
  else if ( strchr (argv[c], 'M') ) PP->problem |= PROBLEM_MAXIMAL;
  else error ("M or C has to be specified", EXIT);
  c++;
  
  while ( argv[c][0] == '-' ){
    switch (argv[c][1]){
      case 'K': II->topk.end = atoi (argv[c+1]);
      break; case 'l': II->lb = atoi (argv[c+1]);
      break; case 'u': II->ub = atoi(argv[c+1]);
      break; case 'U': II->frq_ub = (WEIGHT)atof(argv[c+1]);
      break; case 'w': PP->weight_fname = argv[c+1];
      break; case '#': II->max_solutions = atoi(argv[c+1]);
      break; case ',': II->separator = argv[c+1][0];
      break; case 'Q': PP->outperm_fname = argv[c+1];
      break; default: goto NEXT;
    }
    c += 2;
    if ( argc<c+2 ){ PCE_error (); return; }
  }

  NEXT:;
  PP->SG.fname = argv[c];
  if ( II->topk.end==0 ) PP->dense = II->frq_lb = (WEIGHT)atof(argv[c+1]);
  if ( argc>c+2 ) PP->output_fname = argv[c+2];
}


/******************************************************************/
/******************************************************************/
/******************************************************************/

/* add a vertex to clique and update the degrees of other vertices */
void PCE_add_vertex_to_clique (PROBLEM *PP, QUEUE_INT v, double *sum, double *all){
  QUEUE_INT *x;
  *all += PP->II.itemset.t;
  *sum += PP->occ.list[v];
  QUEUE_ins_ele_ (&PP->II.itemset, v);
  MQUE_FLOOP (PP->SG.edge.v[v], x)
      MALIST_mv (&PP->occ, PP->occ.list[*x] + 1, *x, 0);
}
/* remove a vertex to clique and update the degrees of other vertices */
void PCE_rm_vertex_from_clique (PROBLEM *PP, QUEUE_INT v, double *sum, double *all){
  QUEUE_INT *x;
  *sum -= PP->occ.list[v];
  *all -= PP->II.itemset.t;
  QUEUE_rm_ele_ (&PP->II.itemset, v);
  MQUE_FLOOP (PP->SG.edge.v[v], x)
      MALIST_mv (&PP->occ, PP->occ.list[*x] -1, *x, 0);
}

/* compute the minimum degree in the graph induced by PCE_clq */
QUEUE_INT PCE_min_degree_in_clique (PROBLEM *PP){
  QUEUE_INT *x, m=PP->II.itemset.t+1;
  MQUE_FLOOP (PP->II.itemset, x) ENMIN (m, PP->occ.list[*x]);
  return (m);
}

/* collect minimum degree vertices in PCE_clq */
void PCE_cllect_degree_vertices (PROBLEM *PP, QUEUE_INT d){
  ITEMSET *II = &PP->II;
  QUEUE_INT *x;
  II->add.t = 0;
  MQUE_FLOOP (II->itemset, x)
    if ( PP->occ.list[*x] == d ) QUE_INS (II->add, *x);
}

void PCE_enum_type1_children (PROBLEM *PP, QUEUE_INT min_deg, QUEUE_INT min_add_deg){
  ALIST_ID i, j;
  FLOOP (i, min_add_deg, min_deg)
      MALIST_DO_FORWARD (PP->occ, i, j) QUE_INS (PP->itemcand, j);
}

void PCE_enum_type2_children (PROBLEM *PP, QUEUE_INT min_deg, QUEUE_INT min_add_deg){
  SGRAPH *G = &PP->SG;
  ITEMSET *II = &PP->II;
  QUEUE_INT *x, u, k=0, kk=0;
  QUEUE_ID j;
  
  if ( min_deg < min_add_deg ) return;
  kk=0; MALIST_DO_FORWARD (PP->occ, min_deg, u){
    PP->itemary[u] = 0;
    PP->vecary[kk++] = u;
//    printf ("# %d\n", u);
  }
  qsort_VEC_ID (PP->vecary, kk, 1);

  MQUE_FLOOP (II->itemset, x){
    PP->itemary[*x] -= G->edge.t;  // not to find vertices in II->itemset
    if ( PP->occ.list[*x] == min_deg && k < min_deg){
      QUEUE_BE_LOOP_ (G->edge.v[*x], j, u){
        if ( u<*x ) break;
        PP->itemary[u]++;
//        printf ("#inc++ %d\n", u);
      }
      k++;
    }
  }

  PCE_cllect_degree_vertices (PP, min_deg);
  QUE_INS (II->add, G->edge.t);
  j=0; FLOOP (k, 0, kk){
    u = PP->vecary[k];
    while ( II->add.v[j]<u ) j++;
//    printf ("u=%d, j=%d,%d\n", u, j, II->add.v[j]);
    if ( PP->itemary[u] == j ) QUE_INS (PP->itemcand, u);
  }
}

void PCE_enum_type3_children (PROBLEM *PP, QUEUE_INT min_deg, QUEUE_INT min_add_deg){
  SGRAPH *G = &PP->SG;
  ITEMSET *II = &PP->II;
  ALIST_ID u; 
  int flag=1, vmin=0;
  QUEUE_ID jj=0, j;
  QUEUE_INT *x, *y, k=0, kk=0;

  if ( min_deg+1 < min_add_deg ) return;
  MALIST_DO_FORWARD (PP->occ, min_deg+1, u){
    PP->itemary[u] = 0;
    PP->vecary[kk++] = u;
  }
  qsort_QUEUE_INT (PP->vecary, kk, 1);
  MQUE_FLOOP (II->itemset, x){
    PP->itemary[*x] -= G->edge.t;  // not to find vertices in PCE_clq
    if ( PP->occ.list[*x] == min_deg ){
      jj++;
      if ( flag ){ flag=0; vmin=*x; }
      MQUE_FLOOP (G->edge.v[*x], y){
        if ( *y >= vmin ) break;
        PP->itemary[*y]++;
//        printf ("inc++ %d\n", u);
      }
    }
  }

  MQUE_FLOOP (II->itemset, x){
    if ( PP->occ.list[*x] == min_deg+1 ){
      for ( y=G->edge.v[*x].v+G->edge.v[*x].t-1 ; y>=G->edge.v[*x].v ; y--){
        if ( *y < *x ) break;
        PP->itemary[*y]++;
//        printf ("inc %d\n", u);
      }
      k++;
      if (  k >= min_deg+1 ) break;
    }
  }

  PCE_cllect_degree_vertices (PP, min_deg+1);
  QUE_INS (II->add, G->edge.t);
  j=0; FLOOP (k, 0, kk){
    u = PP->vecary[k];
//    printf ("u=%d, j=%d, jj=%d, kk=%d, vmin=%d\n", u, j, jj, kk, vmin);
    if ( u>vmin ) break;
    while ( II->add.v[j]<u ) j++;
    if ( PP->itemary[u] == j+jj ) QUE_INS (PP->itemcand, u);
  }
}

QUEUE PCE_enum_children (PROBLEM *PP, QUEUE_INT min_deg, QUEUE_INT min_add_deg){
  PP->itemcand.t = 0;
  PCE_enum_type1_children (PP, min_deg, min_add_deg);
  PCE_enum_type2_children (PP, min_deg, min_add_deg);
  PCE_enum_type3_children (PP, min_deg, min_add_deg);
  return ( QUEUE_dup_ (&PP->itemcand) );
}

int PCE_check_maximality (PROBLEM *PP, QUEUE_INT min_add_deg){
  ITEMSET *II = &PP->II;
  QUEUE_INT *x;
  QUEUE_ID i;
  ARY_FILL (PP->itemary, min_add_deg, II->itemset.t+1, 0);
  MQUE_FLOOP (II->itemset, x) PP->itemary[PP->occ.list[*x]]++;
  FLOOP (i, min_add_deg, II->itemset.t+1)
    if ( PP->occ.num[i] > PP->itemary[i] ) return (0);
  return (1);
}

/*************************************************************************/
/* pseudo clique enumeration, main iteration */
/*************************************************************************/
void PCE_iter (PROBLEM *PP, QUEUE_INT v, double sum, double all){
  QUEUE_INT *x, min_deg, min_add_deg;
  QUEUE jump;
  ITEMSET *II = &PP->II;
  
  II->iters++;
  PCE_add_vertex_to_clique (PP, v, &sum, &all);  // add vertex v to the current clique
//  printf ("start:  %d:  ", v);
//  QUEUE_print (&PCE_clq);
  if ( II->itemset.t < 2 ) II->frq = 1.0;
  else II->frq = sum/all;
  min_deg = PCE_min_degree_in_clique (PP);

  min_add_deg = (QUEUE_INT)ceil (II->frq_lb*(all+II->itemset.t) - sum);
  if ( min_add_deg < 0 ) min_add_deg = 0;

//  printf ("%d %d\n", (PCE_problem&1),PCE_check_maximality( min_add_deg ));
  if ( (PP->problem&1) || PCE_check_maximality(PP, min_add_deg ) )
      ITEMSET_output_itemset (II, NULL, 0);
  if ( II->itemset.t == II->ub ) goto END;

  jump = PCE_enum_children (PP, min_deg, min_add_deg); // find all vertices generating children
  MQUE_FLOOP (jump, x) PCE_iter (PP, *x, sum, all); // recursion by adding vertex u
  QUEUE_end (&jump);
  END:;
  PCE_rm_vertex_from_clique (PP, v, &sum, &all);  // remove vertex v from the current clique
}


/****************/
/* main routine */
/****************/
int PCE_main (int argc, char *argv[]){
  QUEUE_INT i;
  PROBLEM PP;
  SGRAPH *G = &PP.SG;

// QUEUE PP->itemcand;   /* queue for candidates */
// PP->itemary:  #vertices of the current clique adjacent to each vertex */
// PP->tmp: used as a temporary array, for itemset
// PP->itemcand: used to keep the list of candidates, temporary (solid one is stored in local variable "jump"

  PROBLEM_init (&PP);
  PCE_read_param (argc, argv, &PP);
if ( ERROR_MES ) return (1);
  G->flag |= LOAD_PERM + LOAD_INCSORT + LOAD_RM_DUP + LOAD_EDGE;
  PP.II.flag |= ITEMSET_ADD;
  PROBLEM_load (&PP);
  PROBLEM_alloc (&PP, G->edge.t, G->edge.t, G->edge.eles, NULL, PROBLEM_ITEMARY +PROBLEM_VECARY +PROBLEM_ITEMCAND);
  MALIST_alloc (&PP.occ, G->edge.t, G->edge.t+1); // element=>
  FLOOP (i, 0, G->edge.t) MALIST_ins_tail (&PP.occ, 0, i, 0);
  
  if ( !ERROR_MES ){
    FLOOP (i, 0, G->edge.t) PCE_iter (&PP, i, 0, 0);
    ITEMSET_last_output (&PP.II);
  }

  PROBLEM_end (&PP);
  return (ERROR_MES?1:0);
}

/*******************************************************************************/
#ifndef _NO_MAIN_
#define _NO_MAIN_
int main (int argc, char *argv[]){
  return (PCE_main (argc, argv));
}
#endif
/*******************************************************************************/

#endif

