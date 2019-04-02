/*
  Image Morphing
  Author: LiangYucheng
*/

#include "Morphing.hpp"
#include <fstream>
#include <sstream>

Morphing::Morphing() {
  cout << "Begin construct!" << endl;
  readSource();
  readTarget();
  loadPoints();
  frames = 11;
  cout << "End of construct!" << endl;
}

Morphing::~Morphing() {}

void Morphing::readSource() {
  cout << "Read Source" << endl;
  source.load_bmp("./data/1.bmp");
}

void Morphing::readTarget() {
  cout << "Read Target" << endl;
  target.load_bmp("./data/2.bmp");
}

void Morphing::loadPoints() {
  // load source points
  ifstream srcPointFile;
  srcPointFile.open("./data/sourcepoints.txt");
  if (srcPointFile.is_open()) {
    string s;
    int x, y;
    while (getline(srcPointFile, s)) {
      stringstream ss(s);
      ss >> x >> y;
      sourcePoints.emplace_back(Point(x, y));
    }
  }
  else {
    cout << "sourcepoints File not exist!" << endl;
  }

  // read target points;
  ifstream tarPointFile;
  tarPointFile.open("./data/targetpoints.txt");
  if (tarPointFile.is_open()) {
    string s;
    int x, y;
    while (getline(tarPointFile, s)) {
      stringstream ss(s);
      ss >> x >> y;
      targetPoints.emplace_back(Point(x, y));
    }
  }
  else {
    cout << "targetpoints File not exist!" << endl;
  }

  // read index
  ifstream indexFile;
  indexFile.open("./data/index.txt");
  if (indexFile.is_open()) {
    string s;
    int x, y, z;
    while (getline(indexFile, s)) {
      vector<int> v;
      stringstream ss(s);
      ss >> x >> y >> z;
      v.emplace_back(x - 1);
      v.emplace_back(y - 1);
      v.emplace_back(z - 1);
      index.emplace_back(v);
    }
  }
  else {
    cout << "index File not exist!" << endl;
  }
}


void Morphing::generateTrianglesFromPoints() {
  cout << "Begin Generate Triangles" << endl;
  // from keypoints to triangles
  for (int i = 0; i < index.size(); i++) {
    Triangle t1(sourcePoints[index[i][0]], sourcePoints[index[i][1]], sourcePoints[index[i][2]]);
    sourceTriangles.emplace_back(t1);

    Triangle t2(targetPoints[index[i][0]], targetPoints[index[i][1]], targetPoints[index[i][2]]);
    targetTriangles.emplace_back(t2);
  }

  // generate mid-process points
  for (int i = 0; i < frames; i++) {
    vector<Point> v;
    for (int j = 0; j < sourcePoints.size(); j++) {
      // Interpolation
      float x = float(sourcePoints[j].x) + float(i + 1) / (frames + 1) * (targetPoints[j].x - sourcePoints[j].x);
      float y = float(sourcePoints[j].y) + float(i + 1) / (frames + 1) * (targetPoints[j].y - sourcePoints[j].y);
      v.emplace_back(Point(x, y));
    }
    midPoints.emplace_back(v);
  }

  // generate mid-process triangles
  for (int i = 0; i < frames; i++) {
    vector<Triangle> v;
    for (int j = 0; j < index.size(); j++) {
      Triangle t(midPoints[i][index[j][0]], midPoints[i][index[j][1]], midPoints[i][index[j][2]]);
      v.emplace_back(t);
    }
    midTriangles.emplace_back(v);
  }
  cout << "End of Generate Triangles" << endl;
}

CImg<float> Morphing::generateMatrix(Triangle oldT, Triangle newT) {
  CImg<float> A(3, 3, 1, 1, 1);
  CImg<float> y1(1, 3, 1, 1, 0);
  CImg<float> y2(1, 3, 1, 1, 0);
  CImg<float> c1(1, 3, 1, 1, 0);
  CImg<float> c2(1, 3, 1, 1, 0);

  // A transposition
  A(0, 0) = oldT.point1.x;
  A(0, 1) = oldT.point2.x;
  A(0, 2) = oldT.point3.x;
  A(1, 0) = oldT.point1.y;
  A(1, 1) = oldT.point2.y;
  A(1, 2) = oldT.point3.y;

  // A' transposition
  y1(0, 0) = newT.point1.x;
  y1(0, 1) = newT.point2.x;
  y1(0, 2) = newT.point3.x;
  y2(0, 0) = newT.point1.y;
  y2(0, 1) = newT.point2.y;
  y2(0, 2) = newT.point3.y;

  c1 = y1.solve(A);
  c2 = y2.solve(A);

  // construct M
  CImg<float> M(3, 3, 1, 1, 0);
  for (int i = 0; i < 3; i++) {
    M(i, 0) = c1(0, i);
    M(i, 1) = c2(0, i);
  }
  M(2, 2) = 1;
  return M;
}

bool Morphing::inTriangle(Point p, Triangle t) {
  float x0 = t.point3.x - t.point1.x;
  float y0 = t.point3.y - t.point1.y;
  float x1 = t.point2.x - t.point1.x;
  float y1 = t.point2.y - t.point1.y;
  float x2 = p.x - t.point1.x;
  float y2 = p.y - t.point1.y;

  float dot00 = x0 * x0 + y0 * y0;
  float dot01 = x0 * x1 + y0 * y1;
  float dot02 = x0 * x2 + y0 * y2;
  float dot11 = x1 * x1 + y1 * y1;
  float dot12 = x1 * x2 + y1 * y2;

  float inverDeno = 1.0 / (dot00 * dot11 - dot01 * dot01);
  float u = (float)(dot11 * dot02 - dot01 * dot12) * inverDeno;
  if (u < 0 || u > 1) {
    return false;
  }

  float v = (dot00 * dot12 -dot01 * dot02) * inverDeno;
  if (v < 0 || v > 1) {
    return false;
  }

  return  u + v <= 1;
}

void Morphing::morphing() {
  cout << "Begin Morphing" << endl;
  int size = midTriangles[0].size();

  // Put origin image into list
  result.push_back(source);

  for (int i = 0; i < frames; i++) {
    CImg<float> mid(target._width, target._height, 1, 3, 1);
    cout << "In frame: " << i+1 << endl;
    cimg_forXY(mid, x, y) {
      CImg<float> x_(1, 3, 1, 1, 1);
      CImg<float> y1(1, 3, 1, 1, 1);
      CImg<float> y2(1, 3, 1, 1, 1);
      for (int m = 0; m < size; m++) {
        Point p(x, y);
        if (inTriangle(p, midTriangles[i][m])) {
          // transform each point in the triangle using Affix Transformation Matrix
          x_(0, 0) = x;
          x_(0, 1) = y;

          // from mid to source
          CImg<float> trans1 = generateMatrix(midTriangles[i][m], sourceTriangles[m]);
          y1 = trans1 * x_;
          // from mid to target
          CImg<float> trans2 = generateMatrix(midTriangles[i][m], targetTriangles[m]);

          y2 = trans2 * x_;

          // linear interpolation
          float a = float(i + 1) / (frames + 1);
          mid(x, y, 0) = (1 - a) * source(y1(0, 0), y1(0, 1), 0) + a * target(y2(0, 0), y2(0, 1), 0);
          mid(x, y, 1) = (1 - a) * source(y1(0, 0), y1(0, 1), 1) + a * target(y2(0, 0), y2(0, 1), 1);
          mid(x, y, 2) = (1 - a) * source(y1(0, 0), y1(0, 1), 2) + a * target(y2(0, 0), y2(0, 1), 2);
          break;
        }
      }
    }
    result.push_back(mid);
  }
  // Put target into the list
  result.push_back(target);
  cout << "End of Morphing" << endl;

  // save the result list
  const string dir_name = "./result/";
  for (int i = 0; i < result.size(); i++) {
    string s = to_string(i + 1);
    s += ".bmp";
    result[i].save_bmp((dir_name + s).c_str());
  }
}

int main() {
  Morphing m;
  m.generateTrianglesFromPoints();
  m.morphing();
  return 0;
}