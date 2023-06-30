#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"

int main(int argc, char **argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;

    bool burble = false ; 
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;

    LAGRAPH_TRY (LAGraph_Random_Init (msg)) ;
    LAGRAPH_TRY (readproblem (&G, NULL,
        true, true, false, GrB_UINT64, false, argc, argv)) ;
    
    GrB_Vector *all_parents, *all_mappings ;
    GrB_Matrix coarsened ;

    GrB_Matrix A = G->A ;

    LAGRAPH_TRY (LAGraph_Coarsen_Matching (&coarsened, &all_parents, &all_mappings, G, LAGraph_Matching_random, 1, 1, 1, 17, msg)) ;
    LAGRAPH_TRY (LAGraph_Matrix_Print (coarsened, LAGraph_COMPLETE, stdout, msg)) ;
    // LAGRAPH_TRY (LAGraph_Vector_Print (all_parents[0], LAGraph_COMPLETE, stdout, msg)) ;
    // LAGRAPH_TRY (LAGraph_Vector_Print (all_mappings[0], LAGraph_COMPLETE, stdout, msg)) ;
    /*
    char msg[1024] ;
    LAGraph_Init (msg) ;
    LAGraph_Random_Init (msg) ;
    GrB_Matrix test = NULL , test2 = NULL ;
    GRB_TRY (LAGraph_Random_Matrix (&test, GrB_BOOL, 3, 5, 0.5, 42, msg)) ;
    GRB_TRY (LAGraph_Random_Matrix (&test2, GrB_BOOL, 5, 3, 0.2, 93, msg)) ;
    GRB_TRY (GrB_transpose (test2, NULL, NULL, test, NULL)) ;
    return (GrB_SUCCESS) ;

    3 3 4
    1 2 1.7
    1 3 0.5
    2 1 1.7
    3 1 0.5
    */
}