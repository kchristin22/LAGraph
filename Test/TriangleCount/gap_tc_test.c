//------------------------------------------------------------------------------
// LAGraph/Test/TriangeCount/ttest.c: test program for LAGraph_tricount
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M

// Usage:  ttest < matrixmarketfile.mtx
//         ttest matrixmarketfile.mtx
//         ttest matrixmarketfile.grb

#include "LAGraph.h"

#define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LAGRAPH_FREE_ALL    \
    GrB_free (&thunk) ;     \
    GrB_free (&C) ;         \
    GrB_free (&M) ;         \
    GrB_free (&A) ;

char *method_name (int method)
{
    switch (method)
    {
        case 0:  return ("minitri:    nnz (A*E == 2) / 3"  ) ;
        case 1:  return ("Burkhardt:  sum ((A^2) .* A) / 6") ;
        case 2:  return ("Cohen:      sum ((L*U) .* A) / 2") ;
        case 3:  return ("Sandia:     sum ((L*L) .* L)"    ) ;
        case 4:  return ("Sandia2:    sum ((U*U) .* U)"    ) ;
        case 5:  return ("SandiaDot:  sum ((L*U') .* L)"   ) ;
        case 6:  return ("SandiaDot2: sum ((U*L') .* U)"   ) ;
        default: return ("invalid method\n") ; abort ( ) ;
    }
}

void print_method (int method)
{
    printf ("%s\n", method_name (method)) ;
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;

    GrB_Matrix A = NULL, C = NULL, M = NULL ;
    #if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
        && ( GxB_IMPLEMENTATION >= GxB_VERSION (3,0,1) )
    GxB_Scalar thunk = NULL ;
    #else
    GrB_Vector thunk = NULL ;     // unused, for LAGRAPH_FREE_ALL
    #endif
    LAGRAPH_OK (LAGraph_init ( )) ;

    int ntrials = 3 ;
    printf ("# of trials: %d\n", ntrials) ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max = LAGraph_get_nthreads ( ) ;
    if (Nthreads [1] == 0)
    {
        // create thread list automatically
        Nthreads [1] = nthreads_max ;
        for (int t = 2 ; t <= nt ; t++)
        {
            Nthreads [t] = Nthreads [t-1] / 2 ;
            if (Nthreads [t] == 0) nt = t-1 ;
        }
    }
    printf ("threads to test: ") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_tic (tic) ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    
    if (argc > 1)
    {
        // Usage:
        //      ./ttest matrixfile.mtx
        //      ./ttest matrixfile.grb

        // read in the file in Matrix Market format from the input file
        char *filename = argv [1] ;
        printf ("matrix: %s\n", filename) ;

        // find the filename extension
        size_t len = strlen (filename) ;
        char *ext = NULL ;
        for (int k = len-1 ; k >= 0 ; k--)
        {
            if (filename [k] == '.')
            {
                ext = filename + k ;
                printf ("[%s]\n", ext) ;
                break ;
            }
        }
        bool is_binary = (ext != NULL && strncmp (ext, ".grb", 4) == 0) ;

        if (is_binary)
        {
            printf ("Reading binary file: %s\n", filename) ;
            LAGRAPH_OK (LAGraph_binread (&C, filename)) ;
        }
        else
        {
            printf ("Reading Matrix Market file: %s\n", filename) ;
            FILE *f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGRAPH_OK (LAGraph_mmread(&C, f));
            fclose (f) ;
        }

    }
    else
    {

        // Usage:  ./ttest < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK (LAGraph_mmread(&C, stdin));
    }

    double t_read = LAGraph_toc (tic) ;
    printf ("\nread A time:     %14.6f sec\n", t_read) ;

    LAGraph_tic (tic) ;
    GrB_Index n, nvals ;
    LAGr_Matrix_nrows (&n, C) ;

    // A = spones (C), and typecast to bool
    LAGr_Matrix_new (&A, GrB_BOOL, n, n) ;
    LAGr_apply (A, NULL, NULL, GxB_ONE_BOOL, C, NULL) ;
    GrB_free (&C) ;

    // M = diagonal mask matrix
    LAGr_Matrix_new (&M, GrB_BOOL, n, n) ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        // M(i,i) = true ;
        LAGr_Matrix_setElement (M, (bool) true, i, i) ;
    }

    // make A symmetric (A = spones (A|A')) and remove self edges (via M)
    LAGr_eWiseAdd (A, M, NULL, GrB_LOR, A, A, LAGraph_desc_otcr) ;
    GrB_free (&M) ;
    LAGr_Matrix_nvals (&nvals, A) ;

    double t_process = LAGraph_toc (tic) ;
    printf ("process A time:  %14.6f sec\n", t_process) ;
    printf ("# of nodes: %lu   number of entries: %lu\n", n, nvals) ;

    //--------------------------------------------------------------------------
    // triangle counting
    //--------------------------------------------------------------------------

    // warmup for more accurate timing, and also print # of triangles
    int64_t ntriangles ;
    LAGraph_tic (tic) ;
    LAGRAPH_OK (LAGraph_tricount (&ntriangles, 5, A)) ;
    printf ("# of triangles: %" PRId64 "\n", ntriangles) ;
    double ttot = LAGraph_toc (tic) ;
    printf ("nthreads: %3d time: %12.6f rate: %6.2f (SandiaDot, one trial)\n",
            nthreads_max, ttot, 1e-6 * nvals / ttot) ;

    double t_best = INFINITY ;
    int method_best = -1 ;
    int nthreads_best = -1 ;


    // kron: input graph: nodes: 134217726 edges: 4223264644
    // fails on methods 3 and 4.

    // try all methods 3 to 6
    // for (int method = 3 ; method <= 6 ; method++)

    // just try all methods 5 and 6
    for (int method = 3 ; method <= 6 ; method++)
    {

        printf ("\nMethod: ") ;
        print_method (method) ;
        if (n == 134217726 && method < 5)
        {
            printf ("kron fails on method %d; skipped\n", method) ;
            continue ;
        }

        for (int t = 1 ; t <= nt ; t++)
        {
            int nthreads = Nthreads [t] ;
            if (nthreads > nthreads_max) continue ;
            LAGraph_set_nthreads (nthreads) ;
            int64_t nt2 ;
            LAGraph_tic (tic) ;
            for (int trial = 0 ; trial < ntrials ; trial++)
            {
                LAGRAPH_OK (LAGraph_tricount (&nt2, method, A)) ;
            }
            double ttot = LAGraph_toc (tic) / ntrials ;

            printf ("nthreads: %3d time: %12.6f rate: %6.2f\n", nthreads, ttot,
                1e-6 * nvals / ttot) ;
            if (nt2 != ntriangles) { printf ("Test failure!\n") ; abort ( ) ; }

            LAGr_log (matrix_name, method_name (method), nthreads, ttot) ;

            if (ttot < t_best)
            {
                t_best = ttot ;
                method_best = method ;
                nthreads_best = nthreads ;
            }
        }
    }

    printf ("\nBest method:\n") ;
    printf ("nthreads: %3d time: %12.6f rate: %6.2f ", nthreads_best, t_best,
        1e-6 * nvals / t_best) ;
    print_method (method_best) ;
    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
}

