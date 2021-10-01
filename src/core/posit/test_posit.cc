#include <iostream>
#include <vector>
#include <singa/core/posit.h>
#include "../../test/gtest/gtest.h"

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


TEST(PositTest, Conversion) {
    
    double d = -291713143525;
    posit_t p = double_to_posit(d);
    int err;
    
    uint32_t p_bit_expected = encode_posit(p.ps, p.es, d, &err);
    if (err) cout << "error: " << err << endl;
    uint32_t p_bit_actual;
    pack_posit(p, &p_bit_actual, &err);
    if (err) cout << "error: " << err << endl;

    EXPECT_EQ(p_bit_actual, p_bit_expected);
    // double d2 = -32832;
    // posit_t p2 = double_to_posit(d2);
    // EXPECT_DOUBLE_EQ(d2, posit_to_double(p2));
    // EXPECT_DOUBLE_EQ(posit_to_double(p), posit_to_double(p2)); // somehow this returns true
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

TEST(PositTest, AddSimple) {
    
    std::vector<double> items = { 0, 1, 2, 4, 8, -1, -2, -4, 8, -8, -0.125, 0.125};
    for (uint i = 0; i < items.size(); i++) {
        for (uint j = 0; j < items.size(); j++) {
            cout << items[i] << " " << items[j] << " " << (items[i] + items[j]) << endl;
            double d = items[i] + items[j];
            int err;
            posit_t p = posit_add(double_to_posit(items[i]), double_to_posit(items[j]));
            uint32_t actual;
            uint32_t expected = encode_posit(DEFAULT_PS, DEFAULT_ES, d, &err);
            pack_posit(p, &actual, &err);
            EXPECT_EQ(actual, expected);
            cout << endl;
        }
    }
}

TEST(PositTest, AddFraction) {
    std::vector<double> items { 1, 1000000000, 0.000000001, -1, -1000000000, -0.000000001, 0.29, 3.14, -3, -12.3, -0.569 };
    for (uint i = 0; i < items.size(); i++) {
        for (uint j = 0; j < items.size(); j++) {
            cout << items[i] << " " << items[j] << " " << (items[i] + items[j]) << endl;
            double d = items[i] + items[j];
            int err;
            posit_t p = posit_add(double_to_posit(items[i]), double_to_posit(items[j]));
            uint32_t actual;
            uint32_t expected = encode_posit(DEFAULT_PS, DEFAULT_ES, d, &err);
            pack_posit(p, &actual, &err);
            EXPECT_EQ(actual, expected);
            cout << endl;
        }
    }
}

TEST(PositTest, AddExponent) {
    // { 1, 1000000000, 0.000000001, -1, -1000000000, -0.000000001 };
    // {, 0.29, 3.14, -3, -12.3, -0.569 }
    int err;
    std::vector<vector<double>> items = { 
        { 1, (unsigned long)1 << 16 }, 
        { 1, (unsigned long)1 << 32 }, 
        { 1, (unsigned long)1 << 48 }, 
        { 1, (unsigned long)1 << 63 }, 
    };
    for (uint i = 0; i < items.size(); i++) {
            cout << items[i][0] << " " << items[i][1] << " " << (items[i][0] + items[i][1]) << endl;
            double d = items[i][0] + items[i][1];
            posit_t p = posit_add(double_to_posit(items[i][0]), double_to_posit(items[i][1]));
            uint32_t actual;
            uint32_t expected = encode_posit(DEFAULT_PS, DEFAULT_ES, d, &err);
            pack_posit(p, &actual, &err);
            EXPECT_EQ(actual, expected);
            cout << endl;
    }
}

TEST(PositTest, Modulo) {
    EXPECT_EQ(0 % 8, 0);
    EXPECT_EQ(0 / 8, 0);
    EXPECT_EQ(-1 % 8, -1);
    EXPECT_EQ(-1 / 8, 0);
    EXPECT_EQ(-32 % 8, 0);
    EXPECT_EQ(-32 / 8, -4);
    EXPECT_EQ(-31 % 8, -7);
    EXPECT_EQ(-31 / 8, -3);
    EXPECT_EQ(-33 % 8, -1);
    EXPECT_EQ(-33 / 8, -4);
    EXPECT_EQ(32 % 8, 0);
    EXPECT_EQ(32 / 8, 4);
    EXPECT_EQ(31 % 8, 7);
    EXPECT_EQ(31 / 8, 3);
    EXPECT_EQ(33 % 8, 1);
    EXPECT_EQ(33 / 8, 4);
}


GTEST_API_ int main(int argc, char **argv) {
    posit_init();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

