#include"gtest/gtest.h"
#include"bdmc_hash.hpp"

TEST(hash, Creation) {
  hash H;
}
TEST(hash, Add) {
  hash H;
  EXPECT_TRUE(H.add_momentum(0.2));
}
TEST(hash, Retrieve) {
  hash H;
  H.add_momentum(0.2);
  EXPECT_TRUE(H.has_momentum(0.2));
  EXPECT_FALSE(H.has_momentum(0.3));
}
TEST(hash, MultiAdd) {
  hash H;
  EXPECT_TRUE(H.add_momentum(0.2));
  EXPECT_FALSE(H.add_momentum(0.2));
}
TEST(hash, OutOfRange) {
  hash H;
  EXPECT_THROW(H.add_momentum(1.), std::range_error);
  EXPECT_THROW(H.add_momentum(1.2), std::range_error);
  EXPECT_THROW(H.add_momentum(-0.1), std::range_error);
  EXPECT_THROW(H.add_momentum(-1.e-20), std::range_error);

  H.add_momentum(1.-std::numeric_limits<double>::epsilon());
  H.add_momentum(0.);
}
TEST(hash, NearlyEquivalentIsEquivalent) {
  hash H;
  EXPECT_TRUE(H.add_momentum(0.5));
  EXPECT_FALSE(H.add_momentum(0.5+std::numeric_limits<float>::epsilon()/10.));
  EXPECT_FALSE(H.add_momentum(0.5-std::numeric_limits<float>::epsilon()/10.));
  EXPECT_TRUE(H.add_momentum(0.5+10*std::numeric_limits<float>::epsilon()));
  EXPECT_TRUE(H.add_momentum(0.5-10*std::numeric_limits<float>::epsilon()));
}
