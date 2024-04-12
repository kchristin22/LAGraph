//-------------------------------------------------------------------------------
// LAGraph/src/test/test_KatzCentrality.c: test cases for LAGraph_KatzCentrality
//-------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

//------------------------------------------------------------------------------


#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;

#define LEN 512
char filename [LEN+1] ;

void test_KatzCentrality (void)
{

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    LAGraph_Init (msg) ;
    GrB_Vector Y = NULL;
    GrB_Matrix A = NULL;

    //--------------------------------------------------------------------------
    // test with the west0067 matrix
    //--------------------------------------------------------------------------

    // create the graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "west0067.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A

    GrB_Vector beta = NULL;
    GrB_Index nrows;
    OK (GrB_Matrix_nrows (&nrows, G->A)) ;
    OK (GrB_Vector_new (&beta, GrB_FP64, nrows)) ;
    OK (GrB_Vector_setElement_FP64 (beta, 1.0, 0)) ;

    // test the algorithm
    OK (LAGraph_KatzCentrality (&Y, G, 0.1, &beta, 1000, (double)(1e-06), NULL, true, msg)) ;

    // print the result
    printf ("\nOutput of LAGraph_KatzCentrality:\n") ;
    OK (LAGraph_Vector_Print (Y, LAGraph_COMPLETE, stdout, msg)) ;

    // check the result (ensure Y is equal to G->A)
    // bool ok ;
    // OK (LAGraph_Matrix_IsEqual (&ok, Y, G->A, msg)) ;
    // TEST_CHECK (ok) ;

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&Y)) ;
    OK (LAGraph_Delete (&G, msg)) ;

    LAGraph_Finalize (msg) ;
}

//----------------------------------------------------------------------------
// the make program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
{
    {"KatzCentrality", test_KatzCentrality},    // just one test in this example
    {NULL, NULL}
} ;
