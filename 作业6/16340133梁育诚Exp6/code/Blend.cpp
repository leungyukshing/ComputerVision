#include "Blend.h"

// 将两张图按照偏移拼接在一起
CImg<float> blending(const CImg<float>& A, const CImg<float>& B, int offset_x, int offset_y, int min_x, int min_y) {
  // 在当前图的左边拼接
  if (offset_x > 0) {
    // 新图片宽度和高度（B较大）
    int nwidth = B._width + abs(offset_x);
    int nheight = B._height + abs(offset_y);

    CImg<float> result(nwidth, nheight, 1, B.spectrum(), 0);

    // 以A为基础构建新图
    cimg_forXY(A, i, j) {
      // 减少不必要的赋值（大于min_x的都是B）
      if (i > min_x)
        continue;

      // 三个色道都要赋值
      for (int k = 0; k < A.spectrum(); k++)
        result(i, j, 0, k) = A(i, j, 0, k);
    }

    // 按照偏移量将B拼接（A左，B右）
    cimg_forXY(B, x, y) {
      if (x + offset_x < 0 || x + offset_x > result._width || y + offset_y < 0 || y + offset_y > result._height)
        continue;
      // 三个色道都要赋值
      for (int k = 0; k < B.spectrum(); k++)
        result(x + offset_x, y + offset_y, 0, k) = B(x, y, 0, k);
    }

    return result;
  }
  // 在当前图的右边拼接（B左，A右）
  else {
    // 新图片宽度和高度（A较大）
    int nwidth = A._width + abs(offset_x);
    int nheight = A._height + abs(offset_y);

    CImg<float> result(nwidth, nheight, 1, A.spectrum(), 0);

    // 以B为基础构建新图
    cimg_forXY(B, i, j) {
      // 三个色道都要赋值
      for (int k = 0; k < B.spectrum(); k++)
        result(i, j, 0, k) = B(i, j, 0, k);
    }
    // 按偏移量将A拼接
    cimg_forXY(A, i, j) {
      // 小于min_x的都是B
      if (i < min_x)
        continue;
      if (i - offset_x < 0 || i - offset_x > result._width || j - offset_y < 0 || j - offset_y > result._height)
        continue;
      // 三个色道都要赋值
      for (int k = 0; k < A.spectrum(); k++) {
        result(i - offset_x, j - offset_y, 0, k) = A(i, j, 0, k);
      }
    }
    return result;
  }
}