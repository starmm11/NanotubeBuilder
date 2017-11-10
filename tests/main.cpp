#include <iostream>
#include "NanotubeBuilder.h"

#include <gtest/gtest.h>
using namespace std;

TEST(Basic, Dirs) {
    NTBuilder::crystal_dir a = {1,1,0};
    NTBuilder::crystal_dir b = {1,1,0};
    NTBuilder::crystal_dir c = NTBuilder::operator+(a, b);
    EXPECT_EQ(c[0], 2);
    EXPECT_EQ(c[1], 2);
    EXPECT_EQ(c[2], 0);
}