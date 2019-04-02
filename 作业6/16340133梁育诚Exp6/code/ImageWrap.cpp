#include "ImageWrap.h"

const double PI = 3.1415926;

// 对图像进行柱面投影预处理
CImg<float> CylinderProjection(CImg<float>& img) {
  CImg<float> result;
  
  result.fill(0.0f);

  int width = img._width;
  int height = img._height;
  int depth = img._depth;
  
  result.assign(width, height, depth, 3);

  float centerX = width / 2;
  float centerY = height / 2;
  float f = width / (2 * tan(PI / 4 / 2));

  cimg_forXY(img, i, j) {
    // 计算曲面投影后的新坐标
    float theta = asin((i - centerX) / f);
    int pointX = (f * tan((i - centerX) / f) + centerX);
    int pointY = ((j - centerY) / cos(theta) + centerY);


    for (int k = 0; k < depth; k++) {
      if (pointX >= 0 && pointX < width && pointY >= 0 && pointY < height) {
        result(i, j, k, 0) = img(pointX, pointY, k, 0);
        result(i, j, k, 1) = img(pointX, pointY, k, 1);
        result(i, j, k, 2) = img(pointX, pointY, k, 2);
      }
      else {
        result(i, j, k, 0) = 0;
        result(i, j, k, 1) = 0;
        result(i, j, k, 2) = 0;
      }
    }
  }
  return result;
}


int getXAfterWarping(float x, float y, MatrixXf& H) {
  return (H(0, 0) * x + H(1, 0) * y + H(2, 0)) / (H(6, 0) * x + H(7, 0) * y + 1);
}

int getYAfterWarping(float x, float y, MatrixXf& H) {
  return (H(3, 0) * x + H(4, 0) * y + H(5, 0)) / (H(6, 0) * x + H(7, 0) * y + 1);
}