#include <iostream>
#include "LatticeOperations.h"
#include <gtest/gtest.h>

using namespace std;

TEST(Basic, Dirs) {
    using namespace NTBuilder;
    crystal_dir a = {1, 1, 0};
    crystal_dir b = {2, 1, 0};
    crystal_dir c = a+b;
    EXPECT_EQ(c[0], 3);
    EXPECT_EQ(c[1], 2);
    EXPECT_EQ(c[2], 0);
}
