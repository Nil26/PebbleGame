//
// Created by Shang Zhang on 5/7/17.
//

#include "gtest/gtest.h"
#include "site.h"
//#include "SiteRP_cont.h"
#include "SiteRP.h"

TEST(SiteRP_cont, Creation)
{
    SiteRP a;
}

TEST(SiteRP_cont, ReadFile)
{
    SiteRP a;
    site bbb(0.1,0.2);
    a.OneTrialTest(0,1);
}