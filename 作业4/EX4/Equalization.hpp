#include <iostream>
#include <cmath>
#include <string>
#include "CImg.h"

using namespace cimg_library;
using namespace std;

class Equalization {
public:
  Equalization(string source);
  Equalization(string source, string target);
  ~Equalization();
  void equalization();
  void colorTransform();

private:
  // Equalization
  int imageRows;
  int imageColumns;
  CImg<float> image;
  CImg<float> output;

  // Color Transformation
  int originRows;
  int originColumns;
  int targetRows;
  int targetColumns;
  CImg<float> origin;
  CImg<float> target;
  CImg<float> result;
};