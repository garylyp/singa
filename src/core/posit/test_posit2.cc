#include <iostream>
#include <vector>
#include "posit2.h"
#include "../../test/gtest/gtest.h"

using namespace std;
using namespace singa;



TEST(PositTest, Equal) {
    EXPECT_TRUE(posit_is_equal(posit_zero, posit_zero));
}

TEST(PositTest, EqualityBySign) { 
    posit_t a = posit_from_float(-123);
    posit_t b = posit_from_float(123);
    EXPECT_TRUE(posit_is_smaller(a, b));
    EXPECT_TRUE(posit_is_smaller_equal(a, b));
    EXPECT_FALSE(posit_is_equal(a, b));
    EXPECT_FALSE(posit_is_bigger(a, b));
    EXPECT_FALSE(posit_is_bigger_equal(a, b));
}

TEST(PositTest, Conversion) {
    
    double d = -291713143525;
    posit_t p = posit_from_float(d);
    EXPECT_DOUBLE_EQ(d, posit_to_float(p));
    cout << posit_to_float(posit_maximum) << endl;
    cout << posit_to_float(posit_minimum) << endl;
    // double d2 = -32832;
    // posit_t p2 = double_to_posit(d2);
    // EXPECT_DOUBLE_EQ(d2, posit_to_double(p2));
    // EXPECT_DOUBLE_EQ(posit_to_double(p), posit_to_double(p2)); // somehow this returns true
}


GTEST_API_ int main(int argc, char **argv) {
    posit_init();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

