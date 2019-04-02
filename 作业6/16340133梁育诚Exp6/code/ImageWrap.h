#ifndef warping_h
#define warping_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <cmath>
#include "CImg.h"

using namespace std;
using namespace Eigen;
using namespace cimg_library;
// 曲面投影及透视变换

extern "C" {
#include <vl/generic.h>
#include <vl/sift.h>
#include <vl/kdtree.h>
}

CImg<float> CylinderProjection(CImg<float>& img);

int getXAfterWarping(float x, float y, MatrixXf& H);

int getYAfterWarping(float x, float y, MatrixXf& H);


#endif
