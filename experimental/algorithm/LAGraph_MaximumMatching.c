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

#include "LAGraphX.h"
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
#define VERTEX_DEFN                                                            \
    "typedef struct "                                                          \
    "{ "                                                                       \
    "uint64_t parentC; "                                                       \
    "uint64_t rootC; "                                                         \
    "} "                                                                       \
    "vertex; "

void *initFrontier(vertex *z, void *x, uint64_t i, uint64_t j, const void *y)
{
    z->parentC = i;
    z->rootC = i;
}

#define INIT_FRONTIER_DEFN                                                     \
    "void *initFrontier(vertex *z, void *x, uint64_t i, uint64_t j, const "    \
    "void "                                                                    \
    "*y) "                                                                     \
    "{ "                                                                       \
    "z->parentC = i; "                                                         \
    "z->rootC = i; "                                                           \
    "} "

void *minparent_1(vertex *z, vertex *x, vertex *y)
{
    *z = x->parentC < y->parentC ? *x : *y;
}

#define MIN_PARENT_1_DEFN                                                      \
    "void *minparent_1(vertex *z, vertex *x, vertex *y) "                      \
    "{ "                                                                       \
    "*z = x->parentC < y->parentC ? *x : *y; "                                 \
    "} "

void *minparent_2(vertex *z, vertex *x, vertex *y)
{
    *z = x->parentC < y->parentC ? *x : *y;
}

#define MIN_PARENT_2_DEFN                                                      \
    "void *minparent_2(vertex *z, vertex *x, vertex *y) "                      \
    "{ "                                                                       \
    "*z = x->parentC < y->parentC ? *x : *y; "                                 \
    "} "

// FIXME: revise GraphBLAS so we can tell it that the select2nd operator
// does not use the 'x' input.
void *select2nd(vertex *z, bool *x, vertex *y)
{
    z->parentC = y->parentC;
    z->rootC = y->rootC;
}

#define SELECT_2ND_DEFN                                                        \
    "void *select2nd(vertex *z, bool *x, vertex *y) "                          \
    "{ "                                                                       \
    "z->parentC = y->parentC; "                                                \
    "z->rootC = y->rootC;"                                                     \
    "} "

void *select1st(vertex *z, vertex *x, bool *y)
{
    z->parentC = x->parentC;
    z->rootC = x->rootC;
}

#define SELECT_1ST_DEFN                                                        \
    "void *select2nd(vertex *z, vertex *x, bool *y) "                          \
    "{ "                                                                       \
    "z->parentC = x->parentC; "                                                \
    "z->rootC = x->rootC;"                                                     \
    "} "

void *keepParents(uint64_t *z, vertex *x) { *z = x->parentC; }

#define KEEP_PARENTS_DEFN                                                      \
    "void *keepParents(uint64_t *z, vertex *x) "                               \
    "{ "                                                                       \
    "*z = x->parentC; "                                                        \
    "} "

void *keepRoots(uint64_t *z, vertex *x) { *z = x->rootC; }

#define KEEP_ROOTS_DEFN                                                        \
    "void *keepRoots(uint64_t *z, vertex *x) "                                 \
    "{ "                                                                       \
    "*z = x->rootC; "                                                          \
    "} "

void *buildfCTuples(vertex *z, uint64_t *x, uint64_t i, uint64_t j,
                    const void *y)
{
    z->parentC = i;
    z->rootC = *x;
}

#define BUILT_FC_TUPLES_DEFN                                                   \
    "void *buildfCTuples(vertex *z, uint64_t *x, uint64_t i, uint64_t j, "     \
    "const void *y) "                                                          \
    "{ "                                                                       \
    "z->parentC = i; "                                                         \
    "z->rootC = *x; "                                                          \
    "} "

void *vertexTypecast(vertex *z, uint64_t *x)
{
    z->parentC = *x;
    z->rootC = *x;
}

#define VERTEX_TYPECAST_DEFN                                                   \
    "void *vertexTypecast(vertex *z, uint64_t *x) "                            \
    "{ "                                                                       \
    "z->parentC = *x; "                                                        \
    "z->rootC = *x; "                                                          \
    "} "

void *setParentsMates(vertex *z, vertex *x, vertex *y)
{
    z->parentC = y->parentC;
    z->rootC = x->rootC;
};

#define SET_PARENTS_MATES_DEFN                                                 \
    "void *setParentsMates(vertex *z, vertex *x, vertex *y) "                  \
    "{ "                                                                       \
    "z->parentC = y->parentC; "                                                \
    "z->rootC = x->rootC; "                                                    \
    "} "

//------------------------------------------------------------------------------
// invert
//------------------------------------------------------------------------------

// This function "inverts" an input vector by swapping its row indices
// and its values, returning the result in an output vector.
//
// For example, for the indices/values of an input vector (in) with 5 entries
// and length 100:
//
//      indices: 0  3  5 42 99
//      values:  4 98  1  3 12
//
// on output, the out vector will contain:
//
//      indices: 4 98  1  3 12
//      values:  0  3  5 42 99
//
// The output vector will normally be jumbled since the values will not appear
// in any particular order.  The method assumes that the input values are in
// range 0 to n-1 where n = length(out). The values in the input vector
// may be duplicated and this argument of the function must be set accordingly.
// Both the in vector and out vector must have the same type (GrB_UINT64).  The
// lengths of the two vectors need not be the same, so long as the indices
// remain in range.  Results are undefined if these conditions do not hold.
//
// The in and out vectors may be aliased.  If not aliased, the input vector is
// cleared of all entries on output.  If in and out are aliased, then the
// inversion is performed in-place.
//
// In SuiteSparse:GraphBLAS, this method takes O(1) time if the in vector is in
// CSC (sparse, by column) format.  Otherwise it can take O(e) time if e =
// nvals(in), because the unpack below will convert the in vector to CSC and
// then unpack it.

#undef LG_FREE_ALL
#define LG_FREE_ALL                                                            \
    {                                                                          \
        LAGraph_Free((void *)&I, NULL);                                        \
        LAGraph_Free((void *)&X1, NULL);                                       \
        LAGraph_Free((void *)&X2, NULL);                                       \
    }

static inline GrB_Info invert_nondestructive(
    GrB_Vector out, // input/output.  On input, only the size and type are
                    // kept; any entries in the 'out' vector cleared.  It is
                    // then replaced with the inversion of the input vector.
    GrB_Vector in,  // input vector, not modified.  There must be no duplicate
                    // values in the input vector.
    char *msg)
{
    bool jumbled = 1;
    GrB_Index *I = NULL;
    GrB_Index *X1 = NULL;
    GrB_Index *X2 = NULL; // not used
    GrB_Index IBytes = 0, XBytes = 0;
    uint64_t nvals = 0;

    // the input and output vectors cannot be the same vector
    //  GRB_TRY(GxB_print (in,0)) ; GRB_TRY(GxB_print (out,0)) ;
    ASSERT(in != out);

    // All input/output vectors must be of type GrB_UINT64.

    GRB_TRY(
        GxB_Vector_unpack_CSC(in, (GrB_Index **)&I, (void **)&X1, &IBytes,
                              &XBytes, NULL, &nvals, &jumbled,
                              NULL)); // the output and input should have no
                                      // duplicates, so the order doesn't matter
    GRB_TRY(GrB_Vector_clear(out)); // clear the output first as a prerequisite
                                    // of the build method
    GRB_TRY(GrB_Vector_build_UINT64(
        out, X1, I, nvals,
        NULL)); // build does not take ownership of the lists I and X,
                // but only copies them, these lists will be given
                // again to the input
                // the input should have no duplicates in the
                // values list, so dups are not handled
    GRB_TRY(GxB_Vector_pack_CSC(in, (GrB_Index **)&I, (void **)&X1, IBytes,
                                XBytes, NULL, nvals, jumbled, NULL));
    //  GRB_TRY(GxB_print (in,0)) ; GRB_TRY(GxB_print (out,0)) ;
}

static inline GrB_Info
invert(GrB_Vector out, // input/output.  Same as invert_nondescructive above.
       GrB_Vector in,  // input vector, empty on output (unless in == out)
       bool dups,      // flag that indicates if duplicates exist in the input
                       // vector's values
       char *msg)
{
    // The input and output vectors can be the same vector
    // that is, in == out is OK.
    // All input/output vectors must be of type GrB_UINT64.

    // the output vector will normally be returned in a jumbled state
    bool jumbled = dups ? 0 : 1; // if there are duplicates, we want the indices
                                 // to be ordered so we can keep the min child
    GrB_Index *I = NULL;         // unpack allocates space for these lists
    GrB_Index *X1 = NULL;
    GrB_Index *X2 = NULL; // not used
    GrB_Index IBytes = 0, XBytes = 0;
    uint64_t nvals = 0;
    // printf("invert, dups %d\n", dups);
    //  GRB_TRY(GxB_print (in,2)) ;
    //  GRB_TRY(GxB_print (out,2)) ;

    // #if LAGRAPH_SUITESPARSE
    GRB_TRY(GxB_Vector_unpack_CSC(in, (GrB_Index **)&I, (void **)&X1, &IBytes,
                                  &XBytes, NULL, &nvals, &jumbled,
                                  NULL)); // #else
    // vanilla case using extractTuples and build:
    // allocate I and X GrB_extractTuples(I, X, in, ...) GrB_build(out, X, I,
    // ...) free I and X: LG_FREE_ALL;
    // #endif
    if (!dups)
    {
        GRB_TRY(GxB_Vector_pack_CSC(out, (GrB_Index **)&X1, (void **)&I, XBytes,
                                    IBytes, NULL, nvals, true, NULL));
    }
    else
    {
        GRB_TRY(GrB_Vector_clear(out));
        GRB_TRY(GrB_Vector_build_UINT64(out, X1, I, nvals, GrB_FIRST_UINT64));
        // build copies the lists so they need to be freed in LG_FREE_ALL
        LG_FREE_ALL;
    }
    //  GRB_TRY(GxB_print (in,2)) ;
    //  GRB_TRY(GxB_print (out,2)) ;
}

static inline GrB_Info
invert_2(GrB_Vector out, // input/output
         GrB_Vector in1, // input vector, empty on output (unless in1 == out)
         GrB_Vector in2, // input vector, empty on output (unless in2 == out)
         bool dups,      // flag that indicates if duplicates exist in the input
                         // vector's values
         char *msg)
{
    // The input vectors cannot be aliased.  However in1==out or in2==out is
    // OK.  The two input vectors must have the same # of entries.
    // All input/output vectors must be of type GrB_UINT64.
    ASSERT(in1 != in2);

    GrB_Index *I = NULL;
    GrB_Index *X1 = NULL;
    GrB_Index *X2 = NULL;
    GrB_Index IBytes = 0, X1Bytes = 0, X2Bytes = 0;
    uint64_t nvals1 = 0, nvals2 = 0;

    //  GRB_TRY(GxB_print (in1,0)) ;
    //  GRB_TRY(GxB_print (in2,0)) ;
    //  GRB_TRY(GxB_print (out,0)) ;

    GRB_TRY(GxB_Vector_unpack_CSC(in1, (GrB_Index **)&I, (void **)&X1, &IBytes,
                                  &X1Bytes, NULL, &nvals1, NULL, NULL));
    LAGraph_Free((void *)&I, NULL);
    GRB_TRY(GxB_Vector_unpack_CSC(in2, (GrB_Index **)&I, (void **)&X2, &IBytes,
                                  &X2Bytes, NULL, &nvals2, NULL, NULL));
    ASSERT(nvals1 == nvals2);
    if (!dups)
    {
        LAGraph_Free((void *)&I, NULL);
        GRB_TRY(GxB_Vector_pack_CSC(out, (GrB_Index **)&X2, (void **)&X1,
                                    X2Bytes, X1Bytes, NULL, nvals2, true,
                                    NULL)); // the values are not ordered,
                                            // so the indices are jumbled
    }
    else
    {
        GRB_TRY(GrB_Vector_clear(out));
        GRB_TRY(GrB_Vector_build_UINT64(out, X2, X1, nvals2, GrB_FIRST_UINT64));
        LG_FREE_ALL;
    }

    //  GRB_TRY(GxB_print (in1,0)) ;
    //  GRB_TRY(GxB_print (in2,0)) ;
    //  GRB_TRY(GxB_print (out,0)) ;
}

GrB_Info check_matching(GrB_Matrix A, GrB_Vector mateC, char *msg);
GrB_Info check_matching(GrB_Matrix A, GrB_Vector mateC, char *msg)
{
    uint64_t nrows, ncols;

    GrB_Matrix_ncols(&ncols, A);
    GrB_Matrix_nrows(&nrows, A);
    GrB_Index nmatched = 0;
    GrB_Vector mateC_dup;
    GrB_Vector_dup(&mateC_dup, mateC);

    // invert to check for dups
    GrB_Index Ibytes = 0, Jbytes = 0, Xbytes = 0;
    bool jumbled = 1;
    GrB_Index *J = NULL, *X = NULL;
    GxB_Vector_unpack_CSC(mateC_dup, (GrB_Index **)&J, (void **)&X, &Jbytes,
                          &Xbytes, NULL, &nmatched, &jumbled, NULL);

    // pack matched values in a matrix
    GrB_Matrix M = NULL;
    bool *val;
    LAGraph_Malloc((void **)&val, nmatched, sizeof(bool), msg);
    for (uint64_t i = 0; i < nmatched; i++)
    {
        val[i] = 1;
    }
    GrB_Matrix_new(&M, GrB_BOOL, nrows, ncols);

    GrB_Matrix_build_BOOL(M, X, J, val, nmatched, NULL);
    LAGraph_Free((void **)&val, msg);
    LAGraph_Free((void **)&X, msg);
    LAGraph_Free((void **)&J, msg);

    // mask with matrix A to check if all edges are present in A
    GrB_Matrix_assign(M, M, NULL, A, GrB_ALL, nrows, GrB_ALL, ncols,
                      GrB_DESC_S);
    GrB_Index nvalsM = 0;
    GrB_Matrix_nvals(&nvalsM, M);
    // if values have been eliminated then edges do not exist in A
    if (nvalsM != nmatched)
    {
        printf("mateC invalid!\n");
        fflush(stdout);
        abort();
    }
    printf("mateC OK!\n");
    fflush(stdout);

    GrB_Matrix_free(&M);
    GrB_free(&mateC_dup);
}

//------------------------------------------------------------------------------
// LAGraph_MaximumMatching
//------------------------------------------------------------------------------

#undef LG_FREE_WORK
#define LG_FREE_WORK                                                           \
    {                                                                          \
        GrB_free(&pathC);                                                      \
        GrB_free(&parentsR);                                                   \
        GrB_free(&Vertex);                                                     \
        GrB_free(&frontierC);                                                  \
        GrB_free(&frontierR);                                                  \
        GrB_free(&initFrontierOp);                                             \
        GrB_free(&I);                                                          \
        GrB_free(&MinParent_1);                                                \
        GrB_free(&MinParent_Monoid_1);                                         \
        GrB_free(&MinParent_2);                                                \
        GrB_free(&MinParent_Monoid_2);                                         \
        GrB_free(&Select2ndOp);                                                \
        GrB_free(&Select1stOp);                                                \
        GrB_free(&MinParent_2nd_Semiring);                                     \
        GrB_free(&MinParent_1st_Semiring);                                     \
        GrB_free(&getParentsOp);                                               \
        GrB_free(&getRootsOp);                                                 \
        GrB_free(&parentsUpdate);                                              \
        GrB_free(&ufrontierR);                                                 \
        GrB_free(&mateR);                                                      \
        GrB_free(&rootsufR);                                                   \
        GrB_free(&pathUpdate);                                                 \
        GrB_free(&rootufRIndexes);                                             \
        GrB_free(&rootsfR);                                                    \
        GrB_free(&rootfRIndexes);                                              \
        GrB_free(&buildfCTuplesOp);                                            \
        GrB_free(&vertexTypecastOp);                                           \
        GrB_free(&setParentsMatesOp);                                          \
        GrB_free(&vr);                                                         \
        GrB_free(&pathCopy);                                                   \
        GrB_free(&currentMatesR);                                              \
    }

#undef LG_FREE_ALL
#define LG_FREE_ALL                                                            \
    {                                                                          \
        LG_FREE_WORK;                                                          \
        GrB_free(&mateC);                                                      \
    }

int LAGraph_MaximumMatching(
    // output:
    GrB_Vector
        *mateC_handle, // mateC(j) = i : Column j of the C subset is matched
                       // to row i of the R subset (ignored on input)
    // input:
    GrB_Matrix A,  // input adjacency matrix, TODO: this should be a LAGraph
                   // of a BIPARTITE kind
    GrB_Matrix AT, // transpose of the input adjacency matrix, necessary to
                   // perform push-pull optimization
                   // if NULL, the push-pull optimization is not performed
    GrB_Vector mateC_init, // input only, not modified, ignored if NULL
    char *msg)
{
    // printf("starting LAGraph_MaximumMatching\n");
    fflush(stdout);
    double tot = LAGraph_WallClockTime();

    // GrB_set (GrB_GLOBAL, true, GxB_BURBLE) ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Vector pathC = NULL;    // make this bitmap/sparse,
                                // if dense I would have to give
                                //  all the entries and make the matrix 1-based
    GrB_Vector parentsR = NULL; // parents of row nodes that are reachable
                                // from paths of the initial column frontier
    GrB_Type Vertex = NULL;
    GrB_Vector frontierC = NULL;
    GrB_Vector frontierR = NULL;
    GrB_IndexUnaryOp initFrontierOp = NULL;
    GrB_Vector I = NULL; // dense vector of 1's
    GrB_BinaryOp MinParent_1 = NULL;
    GrB_BinaryOp MinParent_2 = NULL;
    GrB_Monoid MinParent_Monoid_1 = NULL;
    GrB_Monoid MinParent_Monoid_2 = NULL;
    GrB_BinaryOp Select2ndOp = NULL;
    GrB_BinaryOp Select1stOp = NULL;
    GrB_Semiring MinParent_2nd_Semiring = NULL;
    GrB_Semiring MinParent_1st_Semiring = NULL;
    GrB_UnaryOp getParentsOp = NULL;
    GrB_UnaryOp getRootsOp = NULL;
    GrB_Vector parentsUpdate = NULL; // unmatched rows of R frontier
    GrB_Vector ufrontierR = NULL;    // unmatched rows of R frontier
    GrB_Vector mateR = NULL;         // mateR(i) = j : Row i of the R subset is
                                     // matched to column j of the C subset
    GrB_Vector rootsufR = NULL;
    GrB_Vector pathUpdate = NULL;
    GrB_Vector rootufRIndexes = NULL;
    GrB_Vector rootsfR = NULL;
    GrB_Vector rootfRIndexes = NULL;
    GrB_IndexUnaryOp buildfCTuplesOp = NULL;
    GrB_UnaryOp vertexTypecastOp = NULL;
    GrB_BinaryOp setParentsMatesOp = NULL;
    GrB_Vector vr = NULL;
    GrB_Vector pathCopy = NULL;
    GrB_Vector currentMatesR = NULL;

    GrB_Vector mateC = NULL;

    LG_CLEAR_MSG;

    LG_ASSERT_MSG(mateC_handle != NULL, GrB_NULL_POINTER,
                  "mateC handle is NULL");

    LG_ASSERT_MSG(A != NULL || AT != NULL, GrB_NULL_POINTER,
                  "A matrix is NULL");

    (*mateC_handle) = NULL;

    bool do_pushpull = (AT != NULL) && (A != NULL);

    uint64_t ncols = 0;
    uint64_t nrows = 0;

    if (A != NULL)
    {
        GRB_TRY(GrB_Matrix_ncols(&ncols, A));
        GRB_TRY(GrB_Matrix_nrows(&nrows, A));
    }
    else
    {
        GRB_TRY(GrB_Matrix_nrows(&ncols, AT));
        GRB_TRY(GrB_Matrix_ncols(&nrows, AT));
    }

    GRB_TRY(GrB_Vector_new(&pathC, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&parentsR, GrB_UINT64, nrows));

    GRB_TRY(GxB_Type_new(&Vertex, sizeof(vertex), "vertex", VERTEX_DEFN));

    GRB_TRY(GrB_Vector_new(&frontierC, Vertex, ncols));

    GRB_TRY(GrB_Vector_new(&frontierR, Vertex, nrows));

    GRB_TRY(GxB_IndexUnaryOp_new(&initFrontierOp, (void *)initFrontier, Vertex,
                                 GrB_BOOL, GrB_BOOL, "initFrontier",
                                 INIT_FRONTIER_DEFN));

    GRB_TRY(GrB_Vector_new(&I, GrB_BOOL, ncols));
    GRB_TRY(GrB_Vector_assign_BOOL(I, NULL, NULL, 1, GrB_ALL, ncols, NULL));

    GRB_TRY(GxB_BinaryOp_new(&MinParent_1, (void *)minparent_1, Vertex, Vertex,
                             Vertex, "minparent_1", MIN_PARENT_1_DEFN));
    GRB_TRY(GxB_BinaryOp_new(&MinParent_2, (void *)minparent_2, Vertex, Vertex,
                             Vertex, "minparent_2", MIN_PARENT_2_DEFN));
    vertex infinityParent_1 = {GrB_INDEX_MAX + 1, 0};
    vertex infinityParent_2 = {GrB_INDEX_MAX + 1, 0};
    GRB_TRY(GrB_Monoid_new_UDT(&MinParent_Monoid_1, MinParent_1,
                               &infinityParent_1));
    GRB_TRY(GrB_Monoid_new_UDT(&MinParent_Monoid_2, MinParent_2,
                               &infinityParent_2));

    GRB_TRY(GxB_BinaryOp_new(&Select2ndOp, (void *)select2nd, Vertex, GrB_BOOL,
                             Vertex, "select2nd", SELECT_2ND_DEFN));

    GRB_TRY(GrB_Semiring_new(&MinParent_2nd_Semiring, MinParent_Monoid_2,
                             Select2ndOp));

    GRB_TRY(GxB_BinaryOp_new(&Select1stOp, (void *)select1st, Vertex, Vertex,
                             GrB_BOOL, "select1st", SELECT_1ST_DEFN));

    GRB_TRY(GrB_Semiring_new(&MinParent_1st_Semiring, MinParent_Monoid_1,
                             Select1stOp));

    GRB_TRY(GxB_UnaryOp_new(&getParentsOp, (void *)keepParents, GrB_UINT64,
                            Vertex, "keepParents", KEEP_PARENTS_DEFN));

    GRB_TRY(GxB_UnaryOp_new(&getRootsOp, (void *)keepRoots, GrB_UINT64, Vertex,
                            "keepRoots", KEEP_ROOTS_DEFN));

    GRB_TRY(GrB_Vector_new(&parentsUpdate, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&ufrontierR, Vertex, nrows));

    GRB_TRY(GrB_Vector_new(&rootsufR, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&pathUpdate, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&rootufRIndexes, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&rootsfR, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&rootfRIndexes, GrB_UINT64, ncols));

    GRB_TRY(GxB_IndexUnaryOp_new(&buildfCTuplesOp, (void *)buildfCTuples,
                                 Vertex, GrB_UINT64, GrB_BOOL, "buildfCTuples",
                                 BUILT_FC_TUPLES_DEFN));

    GRB_TRY(GxB_UnaryOp_new(&vertexTypecastOp, (void *)vertexTypecast, Vertex,
                            GrB_UINT64, "vertexTypecast",
                            VERTEX_TYPECAST_DEFN));

    GRB_TRY(GxB_BinaryOp_new(&setParentsMatesOp, (void *)setParentsMates,
                             Vertex, Vertex, Vertex, "setParentsMates",
                             SET_PARENTS_MATES_DEFN));

    GRB_TRY(GrB_Vector_new(&vr, GrB_UINT64, nrows));

    GRB_TRY(GrB_Vector_new(&pathCopy, GrB_UINT64, ncols));

    GRB_TRY(GrB_Vector_new(&currentMatesR, GrB_UINT64, nrows));

    uint64_t npath = 0;
    bool y = 0;
    double mxm_time = 0;

    GRB_TRY(GrB_Vector_new(&mateC, GrB_UINT64, ncols));
    GRB_TRY(GrB_Vector_new(&mateR, GrB_UINT64, nrows));

    if (mateC_init != NULL)
    {
        uint64_t nmatched = 0;

        // mateC = (uint64_t) mateC_init
        GRB_TRY(
            GrB_assign(mateC, NULL, NULL, mateC_init, GrB_ALL, ncols, NULL));
        GRB_TRY(GrB_Vector_nvals(&nmatched, mateC));
        if (nmatched)
        {
            // mateR = invert (mateC), but do not clear the input
            LAGRAPH_TRY(invert_nondestructive(mateR, mateC, msg));
        }
        check_matching(A, mateC, msg);
    }

    /* debug
    GxB_Vector_fprint(mateR, "mateR", GxB_COMPLETE, stdout);
    */

    double ttt = LAGraph_WallClockTime();
    int64_t counter = 0;
    do
    {
        uint64_t nmatched = 0;
        GRB_TRY(GrB_Vector_nvals(&nmatched, mateC));
        printf("\niteration %ld, matched %ld time so far %g\n", counter++,
               nmatched, LAGraph_WallClockTime() - ttt);

        GRB_TRY(GrB_Vector_clear(parentsR));
        // for every col j not matched, assign f(j) = VERTEX(j,j)
        GRB_TRY(GrB_Vector_apply_IndexOp_UDT(
            frontierC, mateC, NULL, initFrontierOp, I, &y, GrB_DESC_RSC));

        /* debug
        GrB_Index C[ncols];
        vertex *V = malloc(ncols * sizeof(vertex));
        GrB_Vector_extractTuples_UDT(C, V, &ncols, frontierC);
        for (int k = 0; k < ncols; k++)
        {
            printf("\nfc (%d) = (%ld, %ld)", (int)C[k], V[k].parentC,
        V[k].rootC);
        }
        GxB_Vector_fprint(mateC, "mateC", GxB_COMPLETE, stdout);
        */

        uint64_t nfC = 0;

        do
        {
            //----------------------------------------------------------------------
            // STEPS 1,2: Explore neighbors of column frontier (one step of
            // BFS), keeping only unvisited rows in the frontierR
            //----------------------------------------------------------------------

            printf(".");
            fflush(stdout);

            double t = LAGraph_WallClockTime();
            if (do_pushpull)
            {
                int32_t kind;
                LAGRAPH_TRY(LG_GET_FORMAT_HINT(frontierC, &kind));
                if (kind == LG_BITMAP || kind == LG_FULL)
                {
                    // the frontierC vector is bitmap or full
                    // pull (vector's values are pulled into A)
                    GRB_TRY(GrB_mxv(frontierR, parentsR, NULL,
                                    MinParent_2nd_Semiring, A, frontierC,
                                    GrB_DESC_RSC));
                    printf("mxv\n");
                }
                else
                {
                    // the frontierC vector is sparse or hypersparse
                    // push (vector's values are pushed to A)
                    // GRB_TRY(GrB_Vector_clear(frontierR));
                    GRB_TRY(GrB_vxm(frontierR, NULL, NULL,
                                    MinParent_1st_Semiring, frontierC, AT,
                                    GrB_DESC_R));
                    printf("vxm\n");
                    GRB_TRY(GrB_Vector_assign(frontierR, parentsR, NULL,
                                              frontierR, GrB_ALL, nrows,
                                              GrB_DESC_RSC));
                }
            }
            // else
            // {
            //     if (A != NULL)
            //     {
            //         // Only the pull method can be used if AT is not given
            //         GRB_TRY(GrB_mxv(frontierR, parentsR, NULL,
            //                         MinParent_2nd_Semiring, A, frontierC,
            //                         GrB_DESC_RSC));
            //     }
            //     else
            //     {
            //         // Only the push method can be used if A is not given
            //         GRB_TRY(GrB_vxm(frontierR, parentsR, NULL,
            //                         MinParent_1st_Semiring, frontierC, AT,
            //                         GrB_DESC_RSC));
            //     }
            // }
            t = LAGraph_WallClockTime() - t;
            mxm_time += t;

            //----------------------------------------------------------------------
            // STEPS 3,4: Select univisited, matched and unmatched row vertices
            //----------------------------------------------------------------------
            // set parents of row frontier
            GRB_TRY(GrB_Vector_apply(
                parentsR, frontierR, NULL, getParentsOp, frontierR,
                GrB_DESC_S)); // use input as mask to only update or insert
                              // parents without deleting the ones not
                              // updated

            // select unmatched rows of the R frontier
            GRB_TRY(GrB_Vector_assign(ufrontierR, mateR, NULL, frontierR,
                                      GrB_ALL, nrows, GrB_DESC_RSC));
            // select matched rows of the R frontier
            GRB_TRY(GrB_Vector_assign(frontierR, mateR, NULL, frontierR,
                                      GrB_ALL, nrows, GrB_DESC_RS));

            // keep only mates of rows in frontierR
            GRB_TRY(GrB_Vector_assign(currentMatesR, frontierR, NULL, mateR,
                                      GrB_ALL, nrows, GrB_DESC_RS));

            /* debug
            uint64_t nvals = 0;
            GxB_Vector_fprint(currentMatesR, "currentMatesR", GxB_COMPLETE,
            stdout); GrB_Index *R = (GrB_Index *)malloc(nrows *
            sizeof(GrB_Index)); vertex *VR = (vertex *)malloc(nrows *
            sizeof(vertex)); GrB_Vector_nvals(&nvals, frontierR);
            GrB_Vector_extractTuples_UDT(R, VR, &nrows, frontierR);
            for (int k = 0; k < nrows; k++)
            {
                printf("\nfr (%d) = (%ld, %ld)", (int)R[k], VR[k].parentC,
            VR[k].rootC);
            }
            GxB_Vector_fprint(parentsR, "pr", GxB_COMPLETE, stdout);
            */

            uint64_t nUfR = 0, nfR = 0;
            GRB_TRY(GrB_Vector_nvals(&nUfR, ufrontierR));
            GRB_TRY(GrB_Vector_nvals(&nfR, frontierR));

            if (nUfR)
            {
                //----------------------------------------------------------------------
                // STEP 5: Store endpoints of newly discovered augmenting paths
                //----------------------------------------------------------------------
                // get roots of unmatched row nodes in the R frontier
                GRB_TRY(GrB_Vector_apply(rootsufR, NULL, NULL, getRootsOp,
                                         ufrontierR, NULL));

                // pathUpdate = invert (rootsufR), but need to handle
                // duplicates
                // printf("pathUpdate = invert (rootsufR)\n");
                LAGRAPH_TRY(invert(pathUpdate, rootsufR, true, msg));

                GRB_TRY(GrB_Vector_assign(
                    pathC, pathUpdate, NULL, pathUpdate, GrB_ALL, ncols,
                    GrB_DESC_S)); // update path without deleting the values
                                  // not updated when GrB_ALL is used, ni is
                                  // the number of rows of the vector

                //----------------------------------------------------------------------
                // STEP 6: Prune vertices in trees yielding augmenting paths
                //----------------------------------------------------------------------
                GRB_TRY(GrB_Vector_clear(rootfRIndexes));

                if (nfR)
                {
                    // get roots of row nodes in the current R frontier
                    GRB_TRY(GrB_Vector_apply(rootsfR, NULL, NULL, getRootsOp,
                                             frontierR, NULL));

                    /* debug
                    GxB_Vector_fprint(rootsfR, "rootsfR", GxB_COMPLETE,
                    stdout);
                    */

                    // keep mates and roots of the R frontier (ordered indices)
                    LAGRAPH_TRY(
                        invert_2(rootfRIndexes, currentMatesR, rootsfR, true,
                                 msg)); // rootsfRIndexes(j) = i, where i
                    // is the col mate of the first
                    // row included in the current R
                    // frontier with a col root of j

                    // keep only col roots that are not included in ufR
                    GRB_TRY(GrB_Vector_assign(rootfRIndexes, pathUpdate, NULL,
                                              rootfRIndexes, GrB_ALL, ncols,
                                              GrB_DESC_RSC));

                    //----------------------------------------------------------------------
                    // STEP 7a (ufrontierR not empty): Move values in the
                    // correct positions for the C frontier
                    //----------------------------------------------------------------------
                    // rootfRIndexes = invert (rootfRIndexes), so that
                    // rootfRIndexes(i) = j, where (i,j) = (parentC, rootC) of
                    // the new frontier C
                    // printf("rootfRIndexes = invert (rootfRIndexes)\n");
                    LAGRAPH_TRY(
                        invert(rootfRIndexes, rootfRIndexes, false, msg));
                }

                //----------------------------------------------------------------------
                // STEP 7b (ufrontierR not empty): Build tuple of (parentC,
                // rootC)
                //----------------------------------------------------------------------
                GRB_TRY(GrB_Vector_clear(frontierC));
                GRB_TRY(GrB_Vector_apply_IndexOp_UDT(frontierC, NULL, NULL,
                                                     buildfCTuplesOp,
                                                     rootfRIndexes, &y, NULL));

                /* debug
                GrB_Index C[ncols];
                vertex *V = malloc(ncols * sizeof(vertex));
                GrB_Vector_extractTuples_UDT(C, V, &ncols, frontierC);
                for (int k = 0; k < ncols; k++)
                {
                    printf("\nfc (%d) = (%ld, %ld)", (int)C[k],
                V[k].parentC, V[k].rootC);
                }
                */
            }
            else
            {
                //----------------------------------------------------------------------
                // STEP 7a (ufrontierR is empty): Set parents of the R frontier
                // to their mates
                //----------------------------------------------------------------------
                // typecast mateR to ensure domain match with frontier R and
                // apply op on frontier to set parents to mates
                GRB_TRY(
                    GrB_Vector_apply(frontierR, NULL, setParentsMatesOp,
                                     vertexTypecastOp, currentMatesR,
                                     NULL)); // fR(i) = (column mate of i,
                                             // rootC) add the structural mask

                //----------------------------------------------------------------------
                // STEP 7b (ufrontierR is empty): Move values in the correct
                // positions for the C frontier
                //----------------------------------------------------------------------
                // invert fr and assign to fC
                // (currentMatesR already contains only the rows of fR)
                LAGRAPH_TRY(
                    invert_2(frontierC, frontierR, currentMatesR, false, msg));
            }

            GRB_TRY(GrB_Vector_nvals(&nfC, frontierC));

        } while (nfC);

        //----------------------------------------------------------------------
        // STEP 8: Augment matching by all augmenting paths discovered in
        // this phase
        //----------------------------------------------------------------------

        GRB_TRY(GrB_Vector_nvals(&npath, pathC));
        uint64_t npathCopy = npath;
        while (npath)
        {
            GRB_TRY(GrB_Vector_nvals(&nmatched, mateC));
            printf(" (%ld,%ld)", npath, nmatched);
            fflush(stdout);
            GxB_print(pathC, 2);
            // vr = invert (pathC), leaving pathC empty
            // pathC doesn't have dup values as it stems from an invertion
            // printf("vr = invert (pathC)\n");
            LAGRAPH_TRY(invert(vr, pathC, false, msg));

            /* debug
            GxB_Vector_fprint(vr, "vr", GxB_COMPLETE, stdout);
            GxB_Vector_fprint(parentsR, "parentsR", GxB_COMPLETE, stdout);
            */

            // assign parents of rows to rows
            GRB_TRY(GrB_Vector_assign(
                vr, vr, NULL, parentsR, GrB_ALL, nrows,
                GrB_DESC_S)); // update the values of vr (descriptor needed
                              // to use mask's structure and not values)

            /* debug
            GxB_Vector_fprint(vr, "vr with updated parents", GxB_COMPLETE,
            stdout);
            */

            // // update mateR:  mateR<vr> = vr
            GRB_TRY(GrB_Vector_assign(mateR, vr, NULL, vr, GrB_ALL, nrows,
                                      GrB_DESC_S));

            // pathC = invert (vr), leaving vr empty (vr may have duplicates
            // after parent assignment)
            // printf("pathC = invert (vr)\n");
            // GRB_TRY(GxB_print(vr, 2));
            LAGRAPH_TRY(invert(pathC, vr, false, msg));
            check_matching(A, pathC, msg);

            // GRB_TRY(GxB_print(pathC, 2));

            /* debug
            GxB_Vector_fprint(pathC, "pathC", GxB_COMPLETE, stdout);
            */

            // keep a copy of the previous row matches of the matched cols
            // that will alter mates
            GRB_TRY(GrB_Vector_assign(pathCopy, pathC, NULL, mateC, GrB_ALL,
                                      ncols, GrB_DESC_RS));

            /* debug
            GxB_Vector_fprint(pathCopy, "pathCopy", GxB_COMPLETE, stdout);
            */

            // update mateC
            GRB_TRY(GrB_Vector_assign(mateC, pathC, NULL, pathC, GrB_ALL, ncols,
                                      GrB_DESC_S));

            // At this point, mateR and mateC are in sync.  That is, they
            // are inversions of each other (mateR == invert(mateC) and
            // mateC == invert (mateR) both hold).

            // swap path and pathCopy
            GrB_Vector temp = pathC;
            pathC = pathCopy;
            pathCopy = temp;

            GRB_TRY(GrB_Vector_nvals(&npath, pathC));

            /* debug
            GxB_Vector_fprint(mateC, "mateC", GxB_COMPLETE, stdout);
            */
            // GrB_Vector gunk = NULL ;
            // GrB_Vector_new (&gunk, GrB_UINT64, ncols) ;
            // GrB_assign (gunk, pathC, NULL, mateC, GrB_ALL, ncols, GrB_DESC_S)
            // ; GxB_print (gunk,2) ; GrB_free (&gunk) ;
            check_matching(A, mateC, msg);
        }

        // compute mateR
        // LAGRAPH_TRY(invert_nondestructive(mateR, mateC, msg));

        npath = npathCopy;
    } while (npath); // only in the first and last iteration should this
                     // condition be false

    /* debug
    GxB_Vector_fprint(mateC, "mateC", GxB_COMPLETE, stdout);
    */

    (*mateC_handle) = mateC;
    LG_FREE_WORK;

    // GrB_set (GrB_GLOBAL, false, GxB_BURBLE) ;
    // tot = LAGraph_WallClockTime() - tot;
    // printf("total time %g, mxm time: %g (%g)\n", tot, mxm_time, mxm_time /
    // tot);
    return (GrB_SUCCESS);
}
