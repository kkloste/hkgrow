#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sparfun

#include <boost/test/unit_test.hpp>

#include "sparfun.h"

BOOST_AUTO_TEST_SUITE(Matrix)
BOOST_AUTO_TEST_CASE(rowcol_subset)
{
    {
        sparserow g;
        int rowptr[] = {0,1,2,3,4};
        int cols[] = {2,1,0,3};
        double vals[] = {1,4,3,2};
        g.n = 4;
        g.m = 4;
        g.ai = (int*)&rowptr;
        g.aj = (int*)&cols;
        g.a = (double*)&vals;
        size_t set[] = {2, 3, 0};
        sparserow* g2 = sf_rowcol_subset(&g, &set[0], 3);
        
        BOOST_CHECK( g2->a[0] == 3. );
        BOOST_CHECK( g2->a[1] == 2. );
        BOOST_CHECK( g2->a[2] == 1. );
        
        BOOST_CHECK( g2->ai[0] == 0 );
        BOOST_CHECK( g2->ai[1] == 1 );
        BOOST_CHECK( g2->ai[2] == 2 );
        BOOST_CHECK( g2->ai[3] == 3 );
        
        BOOST_CHECK( g2->aj[0] == 2 );
        BOOST_CHECK( g2->aj[1] == 1 );
        BOOST_CHECK( g2->aj[2] == 0 );
        
        BOOST_CHECK( g2->n == 3 );
        BOOST_CHECK( g2->m == 3 );
        
        
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Graphs)


BOOST_AUTO_TEST_CASE(core_numbers)
{
    {
        sparserow g;
        int rowptr[] = {0,1,2};
        int cols[] = {1,0};
        int core_nums[] = {0,0};
        g.n = 2;
        g.ai = (int*)&rowptr;
        g.aj = (int*)&cols;
        g.m = 2;
        
        sr_graph_core_numbers(&g, (int*)&core_nums);
        BOOST_CHECK(core_nums[0] == 1);
        BOOST_CHECK(core_nums[1] == 1);
    }
    
    {
        sparserow g;
        int rowptr[] = {0,2,4,6};
        int cols[] = {1,2,0,2,1,2};
        int core_nums[] = {0,0,0};
        g.n = 3;
        g.ai = (int*)&rowptr;
        g.aj = (int*)&cols;
        
        sr_graph_core_numbers(&g, (int*)&core_nums);
        BOOST_CHECK(core_nums[0] == 2);
        BOOST_CHECK(core_nums[1] == 2);
        BOOST_CHECK(core_nums[2] == 2);
    }
    
    {
        sparserow g;
        int rowptr[] = {0,2,5,7,8};
        int cols[] = {1,2,0,2,3,1,2,1};
        int core_nums[] = {0,0,0,0};
        g.n = 4;
        g.ai = (int*)&rowptr;
        g.aj = (int*)&cols;
        
        sr_graph_core_numbers(&g, (int*)&core_nums);
        
        BOOST_CHECK(core_nums[0] == 2);
        BOOST_CHECK(core_nums[1] == 2);
        BOOST_CHECK(core_nums[2] == 2);
        BOOST_CHECK(core_nums[3] == 1);
    }
    
    
}

BOOST_AUTO_TEST_SUITE_END()
