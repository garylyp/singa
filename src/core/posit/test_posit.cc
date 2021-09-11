#include <iostream>
#include <singa/core/posit.h>
#include "gtest/gtest.h"

using namespace std;
using namespace singa;


// posit value: NAR
posit_t dummy_posit_nar() {
    posit_t x = posit_nar;
    return x;
}

// posit value: 0
posit_t dummy_posit_00() {
    posit_t x = posit_zero;
    return x;
}


TEST(PositTest, Equal) {
    EXPECT_TRUE(posit_is_equal(posit_zero, posit_zero));
    EXPECT_TRUE(posit_is_equal(posit_nar, posit_nar));
}

TEST(PositTest, ConvertNanFromPosit) { 
    posit_t x_nar = dummy_posit_nar();
    double actual = posit_to_double(x_nar);
    EXPECT_TRUE(isnan(actual));
}

TEST(PositTest, ConvertNanToPosit) { 
    posit_t actual = double_to_posit(NAN);
    posit_t expected = posit_nar;
    EXPECT_TRUE(posit_is_equal(expected, actual));
}

TEST(PositTest, ConvertZeroFromPosit) { 
    posit_t x_zero = dummy_posit_00();
    double actual = posit_to_double(x_zero);
    EXPECT_DOUBLE_EQ(0.0, actual);
}

TEST(PositTest, ConvertZeroToPosit) { 
    posit_t actual = double_to_posit(0.0);
    posit_t expected = posit_zero;
    EXPECT_TRUE(posit_is_equal(expected, actual));
}

TEST(PositTest, EqualityBySign) { 
    posit_t a = double_to_posit(-123);
    posit_t b = double_to_posit(123);
    EXPECT_TRUE(posit_is_less_than(a, b));
    EXPECT_TRUE(posit_is_less_than_equal(a, b));
    EXPECT_FALSE(posit_is_equal(a, b));
    EXPECT_FALSE(posit_is_greater_than(a, b));
    EXPECT_FALSE(posit_is_greater_than_equal(a, b));
}

TEST(PositTest, EqualityByExp) { 
    posit_t a = double_to_posit(-8);
    posit_t b = double_to_posit(-4);
    EXPECT_TRUE(posit_is_less_than(a, b));
    EXPECT_TRUE(posit_is_less_than_equal(a, b));
    posit_t c = double_to_posit(-0.25);
    posit_t d = double_to_posit(-0.125);
    EXPECT_TRUE(posit_is_less_than(c, d));
    EXPECT_TRUE(posit_is_less_than_equal(c, d));
}

TEST(PositTest, EqualityByFrac) { 
    posit_t a = double_to_posit(-8);
    posit_t b = double_to_posit(-4);
    EXPECT_TRUE(posit_is_less_than(a, b));
    EXPECT_TRUE(posit_is_less_than_equal(a, b));
    posit_t c = double_to_posit(-0.25);
    posit_t d = double_to_posit(-0.125);
    EXPECT_TRUE(posit_is_less_than(c, d));
    EXPECT_TRUE(posit_is_less_than_equal(c, d));
}

GTEST_API_ int main(int argc, char **argv) {
    posit_init();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

