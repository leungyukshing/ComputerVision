/*
  Image Morphing
  Author: LiangYucheng
*/

#ifndef MORPHING_H
#define MORPHING_H

#include <iostream>
#include <string>
#include <vector>
#include "CImg.h"


using namespace std;
using namespace cimg_library;

struct Point {
  float x;
  float y;
  Point(float _x, float _y) : x(_x), y(_y) {
  }
};

struct Triangle {
  Point point1, point2, point3;
  Triangle(Point p1, Point p2, Point p3):point1(p1), point2(p2), point3(p3) {
  }
};

class Morphing
{
public:
  Morphing();
  ~Morphing();
  void readSource(); // read source image
  void readTarget(); // read target image
  void loadPoints(); // read key points in both source and target images
  void generateTrianglesFromPoints(); // genereate triangles with key points
  CImg<float> generateMatrix(Triangle oldT, Triangle newT); // compute Affine Transformation Matrix
  bool inTriangle(Point p, Triangle t); // tell a point is whether in the triangle
  void morphing(); 
private:
  CImg<float> source;
  CImg<float> target;

  vector<Point> sourcePoints;
  vector<Point> targetPoints;
  vector<vector<Point>> midPoints; // mid-process key points

  vector<vector<int>> index; // triangles index

  vector<Triangle> sourceTriangles;
  vector<Triangle> targetTriangles;
  vector<vector<Triangle>> midTriangles; // mid-process triangles

  int frames;
  CImgList<float> result; // result list
};
#endif