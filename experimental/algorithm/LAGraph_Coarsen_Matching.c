//------------------------------------------------------------------------------
// LAGraph_Coarsen_Matching: Coarsen an undirected graph using an edge matching
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vidith Madhu, Texas A&M University

//------------------------------------------------------------------------------

/*
This method is used to coarsen an undirected graph. The coarsening is based on a maximal matching,
which is handled by LAGraph_MaximalMatching.
The inputs to this algorithm are as follows in order:
1. an LAGraph_Graph containing the target graph to coarsen
2. the type of matching to perform (random, heavy, or light)
3. whether to retain the size of the graph when coarsening. If 1, then nodes that are eliminated by a coarsening step
are turned into singletons. If 0, the size of the graph is changed and nodes are explicitly relabeled.
4. whether edges that are combined during a coarsening step should have their edge weights summed (for an unweighted graph, this
counts the number of combined edges)
5. How many coarsening steps to perform
6. Random seed used for maximal matching
7. msg for LAGraph error reporting

There are 2 outputs from the function:
1. A GrB_Matrix of the coarsened graph (if the input adjacency matrix is of type GrB_BOOL or GrB_UINT* or GrB_INT*, it will
have type GrB_INT64. Else, it will have the same type as the input matrix).
2. A list of GrB_Vectors (mapping_result) of length nlevels, where if mapping_result[i][u] = v,
then node u in G_{i} maps to node v in G_{i + 1}, where G_0 is the initial graph. Note that this means 
the length of mapping_result[i] is the number of nodes in G_{i}. If preserve_mapping = true, then there is no need 
for such a result, and a NULL pointer is returned.

More specifically, the coarsening step involves a reduction from a graph G to G', where we use a bijection f from
nodes in G to nodes in G'. We can consider f(u) to be the parent of node u. 
For each edge (u, v) in G, we add an edge (f(u), f(v)) to G' iff f(u) != f(v). In our case,
this bijection is given by the maximal matching, where for every matched edge, one of the endpoints of the edge is the
parent (representative) of both endpoints, and any node not part of a matched edge is its own parent.  
*/

#include "LG_internal.h"
#include "LAGraphX.h"

// #define dbg

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK                                \
{                                                   \
    GrB_free(&E) ;                                  \
    GrB_free(&S) ;                                  \
    GrB_free(&matched_edges) ;                      \
    GrB_free(&E_t) ;                                \
    GrB_free(&S_t) ;                                \
    GrB_free(&edge_parent) ;                        \
    GrB_free(&node_parent) ;                        \
    GrB_free(&full) ;                               \
    LAGraph_Free ((void**)(&rows), msg) ;           \
    LAGraph_Free ((void**)(&cols), msg) ;           \
    LAGraph_Free ((void**)(&vals), msg) ;           \
}                                                   \

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    LAGraph_Delete(&G_cpy, msg) ;                   \
    LAGraph_Free((void**) mapping, msg) ;           \
}                                                   \

int LAGraph_Coarsen_Matching
(
    // outputs:
    GrB_Matrix *coarsened,                  // coarsened adjacency
    GrB_Vector **mapping_result,            // array of parent mappings for each level; if preserve_mapping is true, is NULL
    // inputs:
    LAGraph_Graph G,                        // input graph
    LAGraph_Matching_kind matching_type,    // how to perform the coarsening
    int preserve_mapping,                   // preserve original namespace of nodes
    int combine_weights,                    // whether to sum edge weights or just keep the pattern
    GrB_Index nlevels,                      // #of coarsening levels
    uint64_t seed,                          // used for matching
    char *msg
)
{
    LG_CLEAR_MSG ;
    LAGraph_Graph G_cpy = NULL ;            // used for the IncidenceMatrix function
    GrB_Matrix A = NULL ;                   // resulting adjacency matrix (used for output)
    GrB_Matrix E = NULL ;                   // incidence matrix
    GrB_Matrix E_t = NULL ;                 // transpose of incidence
    GrB_Matrix S = NULL ;                   // S matrix (S[i][j] -> node j maps to node i in coarsened graph)
    GrB_Matrix S_t = NULL ;                 // transpose of S
    GrB_Vector matched_edges = NULL ;       // result of maximal matching
    GrB_Vector edge_parent = NULL ;         // points to parent (representative) node for each edge
    GrB_Vector node_parent = NULL ;         // points to parent (representative) node for each node
    GrB_Vector full = NULL ;                // full vector

    GrB_Vector *mapping = NULL ;            // resulting array of parent mappings (used for output)

    // used to build int64 A matrix if needed
    GrB_Index *rows = NULL ;
    GrB_Index *cols = NULL ;
    int64_t *vals = NULL ;
    GrB_Index nvals, nrows ;
    // check properties (no self-loops, undirected)
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        char typename[LAGRAPH_MAX_NAME_LEN] ;
        GrB_Type type ;
        LG_TRY (LAGraph_Matrix_TypeName (typename, G->A, msg)) ;
        LG_TRY (LAGraph_TypeFromName (&type, typename, msg)) ;

        if ((type == GrB_FP32 || type == GrB_FP64) || type == GrB_INT64) {
            // output will keep the same type as input
            GRB_TRY (GrB_Matrix_dup (&A, G->A)) ;
        } else {
            // output will become int64
            // reasoning: want to prevent overflow from combining edges and accomodate negative edge weights
            #ifdef dbg
                printf("Rebuilding A with GrB_INT64\n");
            #endif

            GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;
            GRB_TRY (GrB_Matrix_nrows (&nrows, G->A)) ;

            LG_TRY (LAGraph_Malloc ((void**)(&rows), nvals, sizeof(GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void**)(&cols), nvals, sizeof(GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void**)(&vals), nvals, sizeof(int64_t), msg)) ;
            // extractTuples casts all entries to int64
            GRB_TRY (GrB_Matrix_extractTuples (rows, cols, vals, &nvals, G->A)) ;

            GRB_TRY (GrB_Matrix_new (&A, GrB_INT64, nrows, nrows)) ;
            GRB_TRY (GrB_Matrix_build (A, rows, cols, vals, nvals, GrB_SECOND_INT64)) ;

            LG_TRY (LAGraph_Free ((void**)(&rows), msg)) ;
            LG_TRY (LAGraph_Free ((void**)(&cols), msg)) ;
            LG_TRY (LAGraph_Free ((void**)(&vals), msg)) ;
        }
    }
    else
    {
        // G is not undirected
        LG_ASSERT_MSG (false, -105, "G must be undirected") ;
    }

    LG_ASSERT_MSG (G->nself_edges == 0, -107, "G->nself_edges must be zero") ;
    
    // make new LAGraph_Graph to use for building incidence matrix
    LG_TRY (LAGraph_New (&G_cpy, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    LG_TRY (LAGraph_Cached_NSelfEdges (G_cpy, msg)) ;

    // set A back (LAGraph_New sets A = NULL)
    A = G_cpy->A ;

    GrB_Index num_nodes ;
    GrB_Index num_edges ;

    GRB_TRY (GrB_Matrix_nvals (&num_edges, A)) ;
    GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;
    num_edges /= 2 ; // since undirected

    GRB_TRY (GrB_Matrix_new (&E_t, GrB_FP64, num_edges, num_nodes)) ;

    GRB_TRY (GrB_Matrix_new (&S_t, GrB_BOOL, num_nodes, num_nodes)) ;

    GRB_TRY (GrB_Vector_new (&edge_parent, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&node_parent, GrB_UINT64, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&full, GrB_BOOL, num_nodes)) ;

    GRB_TRY (GrB_assign (full, NULL, NULL, true, GrB_ALL, num_nodes, NULL)) ;

    if (!preserve_mapping) {
        LG_TRY (LAGraph_Malloc ((void**)(&mapping), nlevels, sizeof(GrB_Vector), msg)) ;
    }

    GrB_Index curr_level = 0 ;

    while (nlevels > 0) {
        // get E
        LG_TRY (LAGraph_Incidence_Matrix (&E, G_cpy, msg)) ;
        GRB_TRY (GrB_Matrix_nvals (&num_edges, A)) ;
        num_edges /= 2 ;
        if (!preserve_mapping) {
            GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;
            // resize structures
            // need to resize E_t, edge_parent, node_parent, full
            GRB_TRY (GrB_Vector_resize (node_parent, num_nodes)) ;
            GRB_TRY (GrB_Vector_resize (full, num_nodes)) ;
        }
        GRB_TRY (GrB_Matrix_resize (E_t, num_edges, num_nodes)) ;
        GRB_TRY (GrB_Vector_resize (edge_parent, num_edges)) ;

        GRB_TRY (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;

        // run maximal matching
        LG_TRY (LAGraph_MaximalMatching (&matched_edges, E, matching_type, seed, msg)) ;

/////////////////////
{
        // make edge_parent
        // want to do E_t * full and get the first entry for each edge (mask output with matched_edges)
        GRB_TRY (GrB_mxv (edge_parent, matched_edges, NULL, GxB_MIN_SECONDI_INT64, E_t, full, GrB_DESC_RS)) ;

        #ifdef dbg
            printf("Printing edge_parent for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Vector_Print (edge_parent, LAGraph_COMPLETE, stdout, msg)) ;
            printf("Printing E for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg)) ;
        #endif

        // now, we have edge_parent (each edge points to its parent node)
        // can do E * edge_parent with min_second to get node_parent
        GRB_TRY (GrB_mxv (node_parent, NULL, NULL, GrB_MIN_SECOND_SEMIRING_UINT64, E, edge_parent, NULL)) ;

        // need to do eWiseAdd with [0...(n - 1)] since some nodes might not touch any matched edges, and we want those nodes to point to self
        // TODO: THIS IS BROKEN
        /*
        METHOD 0: broken
        GRB_TRY (GrB_eWiseAdd (node_parent, node_parent, NULL, GxB_SECONDI_INT64, full, node_parent, GrB_DESC_SC)) ;
        printf("Printing node_parent for level (%lld)\n", curr_level) ;
        LG_TRY (LAGraph_Vector_Print (node_parent, LAGraph_COMPLETE, stdout, msg)) ;
        */

        /*
        // METHOD 1: should work
        GRB_TRY (GxB_eWiseUnion (node_parent, node_parent, NULL,
            GxB_SECONDI_INT64, full, 0, node_parent, 0, GrB_DESC_SC)) ;
        */

        #if 0
        // METHOD 5: should work but is very slow
        for (GrB_Index i = 0; i < num_nodes; i++) {
            if (GxB_Vector_isStoredElement (node_parent, i) == GrB_NO_VALUE) {
                GRB_TRY (GrB_Vector_setElement (node_parent, i, i)) ;
            }
        }
        #endif

        // METHOD 2:
        // node_parent<~node_parent,struct> = 0:n-1
        GrB_apply (node_parent, /* mask: */ node_parent, /* accum: */ NULL,
            GrB_ROWINDEX_INT64, full, (int64_t) 0, GrB_DESC_SC) ;

        #if 0
        // METHOD 3: probably slower
        // suppose ramp = 0:n-1
        GrB_apply (ramp, NULL, NULL, GrB_ROWINDEX_INT64, full, 0, NULL) ;
        // node_parent<~node_parent,struct> = ramp
        GrB_assign (node_parent, node_parent, NULL, ramp, GrB_ALL, num_nodes,
            GrB_DESC_SC) ;
        #endif

        #if 0
        // METHOD 4: fails because accum cannot positional
        // node_parent<~node_parent,struct> += 0, where accum is GxB_SECONDI_INT64
        GrB_assign (node_parent, node_parent, GxB_SECONDI_INT64, 0, GrB_ALL, num_nodes,
            GrB_DESC_SC) ;
        #endif

        #ifdef dbg
            printf("Printing node_parent for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Vector_Print (node_parent, LAGraph_COMPLETE, stdout, msg)) ;
            printf("Printing matched edges for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Vector_Print (matched_edges, LAGraph_COMPLETE, stdout, msg)) ;
        #endif
        
        LG_TRY (LAGraph_Parent_to_S (&S, node_parent, preserve_mapping, msg)) ;

        #ifdef dbg
            printf("Printing S for level (%lld)\n", curr_level) ;
            LG_TRY (LAGraph_Matrix_Print (S, LAGraph_COMPLETE, stdout, msg)) ;
        #endif

        GrB_Index S_rows, S_cols ;

        if (!preserve_mapping) {
            // need to resize S_t
            GRB_TRY (GrB_Matrix_nrows (&S_rows, S)) ;
            GRB_TRY (GrB_Matrix_ncols (&S_cols, S)) ;
//          GRB_TRY (GrB_Matrix_resize (S_t, S_cols, S_rows)) ;
        }
        GrB_Matrix_new (S_t, whatever_type, S_cols, S_rows) ;
        GRB_TRY (GrB_transpose (S_t, NULL, NULL, S, NULL)) ;

        GrB_Semiring semiring = combine_weights ? GrB_PLUS_TIMES_SEMIRING_FP64 : LAGraph_any_one_bool ;

        GRB_TRY (GrB_mxm (S, NULL, NULL, semiring, S, A, NULL)) ;
        if (!preserve_mapping) {
            // resize result
            GRB_TRY (GrB_Matrix_resize (A, S_rows, S_rows)) ;
        }
        GRB_TRY (GrB_mxm (A, NULL, NULL, semiring, S, S_t, NULL)) ;
        GrB_free (&S_t) ;

        G_cpy->A = A ;
        // make nself_edges unknown for delete self edges to work
        G_cpy->nself_edges = LAGRAPH_UNKNOWN ;
        // parent nodes for matched edges will form self-edges; need to delete
        LG_TRY (LAGraph_DeleteSelfEdges (G_cpy, msg)) ;
        A = G_cpy->A ;

        // want to free before we reassign what they point to
        GRB_TRY (GrB_free (&S)) ;
}
/////////////////////

        GRB_TRY (GrB_free (&E)) ;
        GRB_TRY (GrB_free (&matched_edges)) ;

        if (!preserve_mapping){
            // record a deep copy of the current node_parent for the current coarsening level
            GRB_TRY (GrB_Vector_dup (mapping + curr_level, node_parent)) ;
            // printf("passed %d 0x%X\n", curr_level, mapping + (curr_level * sizeof(GrB_Vector))) ;
        }

        nlevels-- ;
        curr_level++ ;
    }
    (*coarsened) = A ;
    (*mapping_result) = mapping ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
