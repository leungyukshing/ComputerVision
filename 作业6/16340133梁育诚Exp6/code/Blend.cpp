#include "Blend.h"

// ������ͼ����ƫ��ƴ����һ��
CImg<float> blending(const CImg<float>& A, const CImg<float>& B, int offset_x, int offset_y, int min_x, int min_y) {
  // �ڵ�ǰͼ�����ƴ��
  if (offset_x > 0) {
    // ��ͼƬ��Ⱥ͸߶ȣ�B�ϴ�
    int nwidth = B._width + abs(offset_x);
    int nheight = B._height + abs(offset_y);

    CImg<float> result(nwidth, nheight, 1, B.spectrum(), 0);

    // ��AΪ����������ͼ
    cimg_forXY(A, i, j) {
      // ���ٲ���Ҫ�ĸ�ֵ������min_x�Ķ���B��
      if (i > min_x)
        continue;

      // ����ɫ����Ҫ��ֵ
      for (int k = 0; k < A.spectrum(); k++)
        result(i, j, 0, k) = A(i, j, 0, k);
    }

    // ����ƫ������Bƴ�ӣ�A��B�ң�
    cimg_forXY(B, x, y) {
      if (x + offset_x < 0 || x + offset_x > result._width || y + offset_y < 0 || y + offset_y > result._height)
        continue;
      // ����ɫ����Ҫ��ֵ
      for (int k = 0; k < B.spectrum(); k++)
        result(x + offset_x, y + offset_y, 0, k) = B(x, y, 0, k);
    }

    return result;
  }
  // �ڵ�ǰͼ���ұ�ƴ�ӣ�B��A�ң�
  else {
    // ��ͼƬ��Ⱥ͸߶ȣ�A�ϴ�
    int nwidth = A._width + abs(offset_x);
    int nheight = A._height + abs(offset_y);

    CImg<float> result(nwidth, nheight, 1, A.spectrum(), 0);

    // ��BΪ����������ͼ
    cimg_forXY(B, i, j) {
      // ����ɫ����Ҫ��ֵ
      for (int k = 0; k < B.spectrum(); k++)
        result(i, j, 0, k) = B(i, j, 0, k);
    }
    // ��ƫ������Aƴ��
    cimg_forXY(A, i, j) {
      // С��min_x�Ķ���B
      if (i < min_x)
        continue;
      if (i - offset_x < 0 || i - offset_x > result._width || j - offset_y < 0 || j - offset_y > result._height)
        continue;
      // ����ɫ����Ҫ��ֵ
      for (int k = 0; k < A.spectrum(); k++) {
        result(i - offset_x, j - offset_y, 0, k) = A(i, j, 0, k);
      }
    }
    return result;
  }
}