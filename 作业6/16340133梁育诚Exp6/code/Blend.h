#ifndef blending_h
#define blending_h
// ÕºœÒªÏ∫œ¿‡
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <Eigen/Dense>
#include "CImg.h"

using namespace std;
using namespace Eigen;
using namespace cimg_library;

CImg<float> blending(const CImg<float>& a, const CImg<float>& b,
                           int offset_x, int offset_y, int min_x, int min_y);

#endif