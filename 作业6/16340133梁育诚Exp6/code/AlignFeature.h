#ifndef align_features_h
#define align_features_h
// RANSAC算法及其相关函数
#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <Eigen/Dense>
#include "CImg.h"

#define NUM_OF_PAIR 4
#define CONFIDENCE 0.99
#define INLINER_RATIO 0.5
#define RANSAC_THRESHOLD 4.0

using namespace std;
using namespace cimg_library;
using namespace Eigen;

extern "C" {
#include <vl/generic.h>
#include <vl/sift.h>
#include <vl/kdtree.h>
}

// 关键点对
struct POINT_PAIR {
  VlSiftKeypoint a;
  VlSiftKeypoint b;
  POINT_PAIR(VlSiftKeypoint _x, VlSiftKeypoint _y) {
    a = _x;
    b = _y;
  }
};

int numberOfIteration(float p, float w, int n);

int random(int min, int max);

vector<POINT_PAIR> getPointPairsFromFeatures(map<vector<float>, VlSiftKeypoint>& feature_a,
    map<vector<float>, VlSiftKeypoint>& feature_b);

vector<int> getIndicesOfInlier(vector<POINT_PAIR>& pairs, MatrixXf& H, set<int>& selected_indices);

MatrixXf getHomographyFromPointPairs(vector<POINT_PAIR>& pairs);

MatrixXf leastSquareSolution(vector<POINT_PAIR>& pairs, vector<int>& idxs);

vector<int> RANSAC(const vector<POINT_PAIR>& pairs);

#endif