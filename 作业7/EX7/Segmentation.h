#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include "CImg.h"
#include <iostream>
#include <vector>

using namespace cimg_library;
using namespace std;

#define EDGE 0
#define NOEDGE 255
#define M_PI 3.14159265358979323846
#define PART 60

// Definition of Point
struct Point {
  int x, y, cnt;
  Point(int _x, int _y, int _cnt): x(_x), y(_y), cnt(_cnt) {}

  bool operator<(const Point &p) const {
    return this->y < p.y;
  }
};

// Definition of Line
struct Line {
  double k, b;
  Line(double _k, double _b): k(_k), b(_b) {}
  int index;
  double distance;
  Line(int _index, double _distance) : index(_index), distance(_distance) {}

  bool operator<(const Line& a) const {
    return this->distance < a.distance;
  }
};

class Segmentation {
public:
  Segmentation(CImg<float> &_img);

  float getThreshold();

  // Compute the histogram of the image
  void getHistogram();

  // Do segmentation
  void segment();

  // Detect Edge
  void getEdge();

  // Initialize Hough Space
  void initHoughSpace();

  // Detect Line from houghspace
  void lineDetect();

  // Find local maximum
  int getPartMax(CImg<float>& img, int &x, int &y);

  // draw the outline of A4
  void drawLines();
  
  // Get crosspoint of the rectangle
  void getCrossPoints();

  // Compute Perspective Transformation Matrix
  vector<CImg<float>> computeTransformMatrix(CImg<float> a4);

  // Perspective Transformation
  void warping(CImg<float> srcImg);

  // Get Binary Image
  void binarization(string path);

private:
  int rows;
  int columns;
  CImg<float> img;
  CImg<float> result;
  CImg<float> hough_space;
  float threshold;
  int histogram[256];

  vector<double> sins; // value of all sin
  vector<double> coss; // value of all cos
  double gradient_threshold;
  double vote_threshold;
  double peak_dis;
  int x_max, y_max;
  vector<pair<int, int>> lines; //set of points
  vector<int> lineCount; // lines' votes set
  vector<int> sortedLineCount; // sorted lines set by its votes
  CImg<unsigned char> afterHough; // after using hough transformation
  CImg<float> showEdge; // image to show the lines detected
  vector<Point> intersections; // CrossPoints of the rectangle
  CImg<float> a4;
};
#endif