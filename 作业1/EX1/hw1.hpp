#define PI 3.1415927
#include "CImg.h"
#include <cmath>

using namespace cimg_library;
using namespace std;

const unsigned char blue[] = {0, 0, 255};
const unsigned char yellow[] = {255, 255, 0};

class MyImage
{
public:
  MyImage() {
    // Read
    SrcImg.load_bmp("1.bmp");
    // Display
    SrcImg.display();
  }

  ~MyImage() {
    SrcImg.display();
    SrcImg.save("2.bmp");
  }

  void changeColor() {
    // white to red
    // black to green
    cimg_forXY(SrcImg, x, y) {
      if (SrcImg(x, y, 0) == 255 && SrcImg(x, y, 1) == 255 && SrcImg(x, y, 2) == 255) {
        SrcImg(x, y, 0) = 255;
        SrcImg(x, y, 1) = 0;
        SrcImg(x, y, 2) = 0;
      }
      else if (SrcImg(x, y, 0) == 0 && SrcImg(x, y, 1) == 0 && SrcImg(x, y, 2) == 0) {
        SrcImg(x, y, 0) = 0;
        SrcImg(x, y, 1) = 255;
        SrcImg(x, y, 2) = 0;
      }
    }
  }

  void drawBlueCircle() {
    SrcImg.draw_circle(50, 50, 30, blue);
  }

  void drawBlueCircleByMe() {
    cimg_forXY(SrcImg, x, y) {
      if ((x-50) * (x-50) + (y-50) * (y-50) <= 900) {
        SrcImg(x, y, 0) = 0;
        SrcImg(x, y, 1) = 0;
        SrcImg(x, y, 2) = 255;
      }
    }
  }

  void drawYellowCircleByMe() {
    cimg_forXY(SrcImg, x, y) {
      if ((x-50) * (x-50) + (y-50) * (y-50) <= 9) {
        SrcImg(x, y, 0) = 255;
        SrcImg(x, y, 1) = 255;
        SrcImg(x, y, 2) = 0;
      }
    }
  }

  void drawYellowCircle() {
    SrcImg.draw_circle(50, 50, 3, yellow);
  }

  void drawBlueLine() {
    double radius = (35 / 180.0) * PI;
    SrcImg.draw_line(0, 0, 100 * cos(radius), 100 * sin(radius), blue);
  }

  void drawBlueLineByMe() {
    double radius = (35 / 180.0) * PI;
    cimg_forXY(SrcImg, x, y) {
      if ((y == (int)(tan(radius) * x)) && (x * x + y * y <= 10000)) {
        SrcImg(x, y, 0) = 0;
        SrcImg(x, y, 1) = 0;
        SrcImg(x, y, 2) = 255;
      }
    }
  }

  void Bresenham() {
    double radius = (35 / 180.0) * PI;
    int xEnd = cos(radius) * 100;
    int yEnd = sin(radius) * 100;

    int dx = xEnd - 0; // x的增量
    int dy = yEnd - 0; // y的增量
    int steps;
    float xIncrement, yIncrement, x = 0, y = 0;

    // 谁的增量大
    if (abs(dx) > abs(dy)) {
      steps = abs(dx);
    }
    else {
      steps = abs(dy);
    }

    xIncrement = float(dx) / float(steps); // x每步骤的增量
    yIncrement = float(dy) / float(steps); // y每步骤增量
    for (int k = 0; k <= steps; k++) {
      x += xIncrement;
      y += yIncrement;
      SrcImg.draw_point(x, y, blue);
    }
  }
private:
  CImg<unsigned char> SrcImg;
};

