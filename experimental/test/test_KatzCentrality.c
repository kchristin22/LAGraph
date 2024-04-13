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

double west_katz [] = {
    1.00014139732199,
    0.98113527961647,
    0.990961312789564,
    0.963092330051061,
    0.982766337720939,
    0.981824600944108,
    0.977267535976225,
    0.986355750024111,
    0.992284041040597,
    0.740129885880301,
    0.905941809311526,
    0.909227587602699,
    0.913767603829989,
    0.918157860803583,
    0.982307449781511,
    0.996030816773654,
    0.990907491072482,
    0.988976762414143,
    0.988702860970875,
    0.971222043176167,
    1.04395316623434,
    1.00589622788841,
    0.981511655693586,
    0.957471456340425,
    1.00544419096611,
    1.00086260153766,
    0.999349601338023,
    1.00042649251723,
    0.9908607924208,
    0.62713733298148,
    0.995559950128439,
    0.98190858689575,
    0.996753532822644,
    0.996821668930618,
    0.997673795353486,
    1.07800218546814,
    1.08577403456582,
    1.03271864944444,
    0.985887443789476,
    0.98995587336602,
    1.01013312291478,
    1.00570122286564,
    1.00584055919791,
    0.995049455388604,
    0.394255898349557,
    1.17249749165967,
    1.15547230199532,
    1.024107770888,
    0.986944071096567,
    1.00924102188452,
    0.943588120770824,
    0.976821954441915,
    0.983098648229359,
    0.965912397153784,
    0.329191049957156,
    1.09887050171704,
    1.45339405870544,
    1.48011665699141,
    1.49942818468319,
    1.50511665750309,
    1.48997808881957,
    1.50243911082964,
    1.4878677116125,
    1.46186370003113,
    1.70369241034709,
    1.47413819664049,
    1.76299950664459
};

double west_katz_norm [] = 
{
    0.112580853338478,
    0.110441430897143,
    0.111547497701818,
    0.108410427416779,
    0.110625030849847,
    0.110519024308932,
    0.110006058578109,
    0.111029072819469,
    0.111696390524008,
    0.083312673944731,
    0.101977282652625,
    0.102347146078829,
    0.102858192719241,
    0.103352381718702,
    0.110573376157874,
    0.112118146098194,
    0.111541439263627,
    0.111324106913919,
    0.111293275215225,
    0.109325548062192,
    0.117512522344109,
    0.11322864547839,
    0.110483797646536,
    0.107777714121893,
    0.113177761970728,
    0.112662035645548,
    0.112491724873459,
    0.112612945061157,
    0.111536182632845,
    0.0705936743509696,
    0.11206514302396,
    0.110528478182284,
    0.112199498584687,
    0.112207168321386,
    0.112303087878456,
    0.121345248047638,
    0.122220085751364,
    0.116247909669901,
    0.110976357957691,
    0.111434320476492,
    0.11370557129993,
    0.113206694750295,
    0.113222379135754,
    0.112007679215753,
    0.0443793903429752,
    0.131982106231898,
    0.130065666830636,
    0.115278626668044,
    0.111095297143901,
    0.113605152003699,
    0.106214937328736,
    0.10995590172079,
    0.110662437361278,
    0.108727868092414,
    0.0370553697855131,
    0.123694288750611,
    0.16360121059308,
    0.16660923818449,
    0.168783039081675,
    0.169423361666033,
    0.167719289636541,
    0.169121963793534,
    0.167481735159269,
    0.164554595234962,
    0.191776028766274,
    0.165936273171984,
    0.198451928322155
};

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
    OK (LAGraph_KatzCentrality (&Y, G, 0.1, &beta, 1000, (double)(1e-06), NULL, false, msg)) ;

    // print the result
    printf ("\nOutput of LAGraph_KatzCentrality:\n") ;
    OK (LAGraph_Vector_Print (Y, LAGraph_COMPLETE, stdout, msg)) ;

    // check the result 
    bool ok = true;
    for (int i = 0; i < nrows; i++){
        double value ;
        OK(GrB_Vector_extractElement_FP64 (&value, Y, i));
        if (abs(value - west_katz[i]) > 1e-05)
        { 
            ok = false ;
            break;
        }
    }
    TEST_CHECK (ok) ;

    // test the algorithm
    OK (LAGraph_KatzCentrality (&Y, G, 0.1, &beta, 1000, (double)(1e-06), NULL, true, msg)) ;

    // print the result
    printf ("\nOutput of LAGraph_KatzCentrality, normalized:\n") ;
    OK (LAGraph_Vector_Print (Y, LAGraph_COMPLETE, stdout, msg)) ;

    // check the result 
    ok = true;
    for (int i = 0; i < nrows; i++){
        double value ;
        OK(GrB_Vector_extractElement_FP64 (&value, Y, i));
        if (abs(value - west_katz_norm[i]) > 1e-05)
        { 
            ok = false ;
            break;
        }
    }
    TEST_CHECK (ok) ;

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&Y)) ;
    OK (GrB_free (&beta)) ;
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
