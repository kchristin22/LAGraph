//------------------------------------------------------------------------------
// LAGraph_MaximumMatching: maximum matching between nodes of disjoint sets
//                          in bipartite graphs
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Christina Koutsou, Aristotle University of Thessaloniki

// add paper

//------------------------------------------------------------------------------

// add explanation of paper

// #define LG_FREE_WORK \
//     {
// }

// #define LG_FREE_ALL \
//     {               \
//         // LG_FREE_WORK ;                  \
//     // GrB_free (centrality) ;
// }

#include "LG_internal.h"

//------------------------------------------------------------------------------
// the Vertex tuple: (parentC, rootC)
//------------------------------------------------------------------------------

typedef struct
{
    uint64_t parentC;
    uint64_t rootC;
} vertex;

// repeat the typedef as a string, to give to GraphBLAS
#define VERTEX_DEFN      \
    "typedef struct "    \
    "{ "                 \
    "uint64_t parentC; " \
    "uint64_t rootC; "   \
    "} "                 \
    "vertex; "

int a = 0;
void *initFrontier(vertex *z, void *x, uint64_t i, uint64_t j, const void *y)
{
    z->parentC = i;
    z->rootC = i;
}

#define INIT_FRONTIER_DEFN                                                           \
    "void *initFrontier(vertex *z, void *x, uint64_t i, uint64_t j, const void *y) " \
    "{ "                                                                             \
    "z->parentC = i; "                                                               \
    "z->rootC = i; "                                                                 \
    "} "

void *minparent(vertex *z, vertex *x, vertex *y)
{
    z = x->parentC < y->parentC ? x : y;
}

#define MIN_PARENT_DEFN                                 \
    "void *minparent(vertex *z, vertex *x, vertex *y) " \
    "{ "                                                \
    "z = x->parentC < y->parentC ? x : y; "             \
    "} "

void *select2nd(vertex *z, bool *x, vertex *y)
{
    z->parentC = y->parentC;
    z->rootC = y->rootC;
}

#define SELECT_2ND_DEFN                               \
    "void *select2nd(vertex *z, bool *x, vertex *y) " \
    "{ "                                              \
    "z->parentC = y->parentC; "                       \
    "z->rootC = y->rootC;"                            \
    "} "

void *keepParents(uint64_t *z, vertex *x)
{
    *z = x->parentC;
}

#define KEEP_PARENTS_DEFN                        \
    "void *keepParents(uint64_t *z, vertex *x) " \
    "{ "                                         \
    "*z = x->parentC; "                          \
    "} "

int LAGraph_MaximumMatching(
    // output/input:
    GrB_Vector *mateC, // mateC(j) = i : Column j of C subset is matched to row i of R subset
    // input:
    LAGraph_Graph G, // input graph
    char *msg)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Matrix A = G->A;
    uint64_t ncols = 0;
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));

    uint64_t nrows = 0;
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));

    GrB_Vector pathC = NULL; // make this bitmap, if dense I would have to give all the entries and make the matrix 1-based
    GRB_TRY(GrB_Vector_new(&pathC, GrB_UINT64, ncols));

    GrB_Vector parentsR = NULL; // parents of row nodes that are reachable from paths of the initial column frontier
    GRB_TRY(GrB_Vector_new(&parentsR, GrB_UINT64, nrows));

    GrB_Type Vertex = NULL;
    GRB_TRY(GxB_Type_new(&Vertex, sizeof(vertex), "vertex", VERTEX_DEFN));

    GrB_Vector frontierC = NULL;
    GRB_TRY(GrB_Vector_new(&frontierC, Vertex, ncols));
    GrB_Vector frontierR = NULL;
    GRB_TRY(GrB_Vector_new(&frontierR, Vertex, nrows));

    GrB_IndexUnaryOp initFrontierOp = NULL;
    GRB_TRY(GxB_IndexUnaryOp_new(&initFrontierOp, (void *)initFrontier, Vertex, GrB_BOOL, GrB_BOOL, "initFrontier", INIT_FRONTIER_DEFN));

    uint64_t nvals = 0;
    bool y = 0; // see if I can get rid of this

    GrB_Vector I = NULL; // dense matrix of 1's
    GRB_TRY(GrB_Vector_new(&I, GrB_BOOL, ncols));
    GRB_TRY(GrB_Vector_assign_INT32(I, NULL, NULL, 1, GrB_ALL, ncols, NULL));

    GrB_BinaryOp MinParent = NULL;
    GRB_TRY(GxB_BinaryOp_new(&MinParent, (void *)minparent, Vertex, Vertex, Vertex, "minparent", MIN_PARENT_DEFN));
    vertex infinityParent = {GrB_INDEX_MAX + 1, 0};
    GrB_Monoid AddMonoid;
    GRB_TRY(GrB_Monoid_new_UDT(&AddMonoid, MinParent, &infinityParent));

    GrB_BinaryOp MultMonoid;
    GRB_TRY(GxB_BinaryOp_new(&MultMonoid, (void *)select2nd,
                             Vertex, GrB_BOOL, Vertex, "select2nd", SELECT_2ND_DEFN));
    GrB_Semiring semiring = NULL;
    GRB_TRY(GrB_Semiring_new(&semiring, AddMonoid, MultMonoid));

    GrB_UnaryOp parentsOp = NULL;
    GRB_TRY(GxB_UnaryOp_new(&parentsOp, (void *)keepParents, GrB_UINT64, Vertex, "keepParents", KEEP_PARENTS_DEFN));

    do
    {
        GRB_TRY(GrB_Vector_clear(pathC));
        GRB_TRY(GrB_Vector_clear(parentsR));
        // for every col j not matched, assign f(j) = VERTEX(j,j)
        GRB_TRY(GrB_Vector_apply_IndexOp_UDT(frontierC, *mateC, NULL, initFrontierOp, I, &y, GrB_DESC_RSC));
        GrB_Index C[ncols];
        vertex *X = malloc(ncols * sizeof(vertex));
        GrB_Vector_extractTuples_UDT(C, X, &ncols, frontierC);
        for (int k = 0; k < ncols; k++)
        {
            printf("\nfc (%d) = (%ld, %ld)", (int)C[k], X[k].parentC, X[k].rootC);
        }
        GxB_Vector_fprint(*mateC, "mateC", GxB_COMPLETE, stdout);

        // do
        // {
        GRB_TRY(GrB_mxv(frontierR, parentsR, NULL, semiring, A, frontierC, GrB_DESC_RSC)); // perform one step of BFS from C nodes and keep only unvisited rows
        GRB_TRY(GrB_Vector_apply(parentsR, NULL, NULL, parentsOp, frontierR, NULL));       // set parents of row frontier

        GRB_TRY(GrB_Vector_nvals(&nvals, frontierC));

        // } while (nvals);

        GrB_Index R[nrows];
        vertex *Y = malloc(nrows * sizeof(vertex));
        GrB_Vector_extractTuples_UDT(R, Y, &nrows, frontierR);
        for (int k = 0; k < nrows; k++)
        {
            printf("\nfr (%d) = (%ld, %ld)", (int)R[k], Y[k].parentC, Y[k].rootC);
        }

        uint64_t P[nrows];
        GrB_Index N[nrows];
        GrB_Vector_extractTuples_UINT64(N, P, &nrows, parentsR);
        for (int k = 0; k < nrows; k++)
        {
            printf("\npr (%d) = %ld", (int)N[k], P[k]);
        }

        LG_ASSERT_MSG(1 == 0, GrB_INVALID_VALUE, "dummy");

        GRB_TRY(GrB_Vector_nvals(&nvals, pathC));

    } while (nvals); // only in the first and last iteration should this condition be false

    return (GrB_SUCCESS);
}