#ifndef stitching_h
#define stitching_h
// ÕºœÒ∆¥Ω”¿‡
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <Eigen/Dense>
#include "CImg.h"
#include "AlignFeature.h"

#define RESIZE_SIZE 500.0
#define THRESHOLD 50.0

using namespace std;
using namespace Eigen;
using namespace cimg_library;

extern "C" {
#include <vl/generic.h>
#include <vl/sift.h>
}

map<vector<float>, VlSiftKeypoint> extractFeatures(CImg<float>& img);

CImg<float> get_gray_image(CImg<float>& img);

void ReplacePairs(vector<POINT_PAIR>& bigger, vector<POINT_PAIR>& smaller);

vector<int> getAvgOffset(const vector<POINT_PAIR>& pairs, vector<int>& idxs);

CImg<float> stitching(vector<CImg<float> >& src_imgs);


#endif