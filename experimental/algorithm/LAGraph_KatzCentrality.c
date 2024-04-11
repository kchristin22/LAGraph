//------------------------------------------------------------------------------
// LAGraph_KCoreDecompose: Helper method to LAGraph_KCore and LAGraph_AllKCore
// that performs graph decomposition given a specified value k.
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Christina Koutsou, Aristotle University of Thessaloniki (AUTh)

//------------------------------------------------------------------------------

#define LG_FREE_WORK                        \
{                                           \
    /* free any workspace used here */      \
    /* GrB_free (&b) ; */                     \
    /*GrB_Matrix_free (&A) ; */                         \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (&C) ;                         \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_KatzCentrality
(
    // outputs:
    GrB_Vector *central,    // katz centrality
    // inputs:
    LAGraph_Graph G,        // input graph
    double alpha,           // attenuation factor
    GrB_Vector *beta,       // bias factor
    uint64_t max_iter,      // maximum number of iterations
    double tol,             // tolerance for convergence
    GrB_Vector *nstart,     // starting vector of centrality for each node
    bool normalize,         // if true, normalize the result
    char *msg
)
{
    // clear user message
    LG_CLEAR_MSG ;

    GrB_Vector C = NULL ;
    LG_ASSERT (central != NULL, GrB_NULL_POINTER) ; // check size of central vector as well?
    (*central) = NULL ;

    // add assert for alpha value
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    GrB_Matrix A = NULL ;
    A = G->A ;

    GrB_Type A_type;
    GRB_TRY(GxB_Matrix_type(&A_type, A));
    LG_ASSERT_MSG (A_type == GrB_FP64, GrB_DOMAIN_MISMATCH, "Adjacency matrix must be of type double") ;

    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    GRB_TRY (GrB_Vector_new (&C, GrB_FP64, n)) ;
    if (nstart == NULL)
    {
        GRB_TRY (GrB_Vector_assign_FP64 (C, NULL, NULL, 0, GrB_ALL, n, NULL)) ; // set all values to 0
    }
    else 
    {
        GrB_Index nstart_size ;
        GRB_TRY (GrB_Vector_size (&nstart_size, *nstart)) ;
        LG_ASSERT_MSG (nstart_size == n, GrB_INVALID_VALUE,
            "nstart vector must have the same size as the number of nodes") ;

        GRB_TRY (GrB_Vector_dup (&C, *nstart)) ; // copy nstart to C
    }

    GrB_Index nvals ; // number of non-zero values in beta vector

    GrB_Vector b;

    GRB_TRY (GrB_Vector_nvals (&nvals, *beta)) ;
    if (nvals == 1){
        // check if needed
        // beta is a scalar
        // Try GrB_Scalar b = *beta;

        // get index of beta value (or assume it's in the first position)
        GxB_Iterator beta_it ;
        GRB_TRY (GxB_Iterator_new (&beta_it)) ;
        GRB_TRY (GxB_Vector_Iterator_attach(beta_it, *beta, NULL));
        GrB_Index beta_value_index = GxB_Vector_Iterator_getp(beta_it);

        // get beta value as scalar
        GRB_TRY (GrB_Vector_extractElement_Scalar (b, *beta, beta_value_index)) ;
    }
    else if (nvals == n)
    {
        // beta is a vector
        GRB_TRY (GrB_Vector_dup (&b, *beta)) ; // copy beta to b
    }
    else
    {
        LG_ASSERT_MSG (false, GrB_INVALID_VALUE, "beta vector must be a scalar or have the same size as the number of nodes") ;
    }

    // if matrix is iso-valued, then we can multiply the iso value with alpha and use plus_first (choose x)
    bool iso_valued = false ;
    GRB_TRY (GxB_Matrix_iso(&iso_valued, A)) ;

    GrB_Semiring semiring ;
    if (iso_valued){
        void *iso_val ;
        GRB_TRY(GxB_Matrix_pack_CSR(A, NULL, NULL, &iso_val, 0, 0, 8, iso_valued, NULL, NULL)) ; // A is of type double

        alpha *=  *(double *)iso_val;

        semiring = GxB_PLUS_SECOND_FP64 ;
    }
    else
    {
        semiring = GxB_PLUS_TIMES_FP64 ;
    }

    GrB_Vector C_prev = NULL ;
    GRB_TRY (GrB_Vector_new (&C_prev, GrB_FP64, n)) ;
    for (uint64_t i = 0; i < max_iter; i++)
    {
        // C_prev = C
        GRB_TRY (GrB_Vector_dup (&C_prev, C)) ;

        // C = alpha * semiring(A,C_prev) + b
        GRB_TRY (GrB_mxv (C, NULL, NULL, semiring, A, C_prev, NULL)) ;
        GRB_TRY (GrB_Vector_apply_BinaryOp1st_FP64 (C, NULL, NULL,  GrB_TIMES_FP64, alpha, C, NULL)) ;
        GRB_TRY (GrB_Vector_eWiseAdd_BinaryOp (C, NULL, NULL, GrB_PLUS_FP64, C, b, NULL)) ;


        // check for convergence
        GrB_Vector diff = NULL ;
        GRB_TRY (GrB_Vector_new (&diff, GrB_FP64, n)) ;
        GRB_TRY (GrB_eWiseAdd (diff, NULL, NULL, GrB_MINUS_FP64, C, C_prev, NULL)) ;
        GRB_TRY (GrB_Vector_apply (diff, NULL, NULL, GrB_ABS_FP64, diff, NULL)) ;
        double norm_diff ;
        GRB_TRY (GrB_Vector_reduce_FP64 (&norm_diff, NULL, GrB_PLUS_MONOID_FP64, diff, NULL)) ;
        if (norm_diff < n * tol)
        {
            if (normalize)
            {
                // normalize C
                double euclidean_norm ;
                GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT8 (C, NULL, NULL, GxB_POW_FP64, C, 2, NULL)) ;
                GRB_TRY (GrB_Vector_reduce_FP64 (&euclidean_norm, NULL, GrB_PLUS_MONOID_FP64, C, NULL)) ;
                GRB_TRY (GrB_Vector_apply_BinaryOp2nd_FP64 (C, NULL, NULL, GrB_DIV_FP64, C, euclidean_norm, NULL)) ;
            }
            (*central) = C ;
            LG_FREE_WORK ;
            return (GrB_SUCCESS) ;
        }
    }

    LG_FREE_ALL ; 
    LG_ASSERT_MSG (false, GxB_EXHAUSTED, "Katz centrality did not converge") ;
}