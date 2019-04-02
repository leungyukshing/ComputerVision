#include <iostream>
#include "Segmentation.h"

using namespace std;

int main() {
  cout << "Image Segmentation" << endl;
  string root = "set1/";

  for (int i = 1; i <= 6; i++) {
    CImg<float> img((root + to_string(i) + ".bmp").c_str());
    CImg<float> grayImg(img._width, img._height, 1, 1);
    
    // Change to Grey Image
    cimg_forXY(img, x, y) {
      grayImg(x, y, 0, 0) = 0.299 * img(x, y, 0, 0) + 0.587 * img(x, y, 0, 1) + 0.114 * img(x, y, 0, 2);
    }

    Segmentation seg(grayImg);
    seg.segment();
    seg.binarization("Seg" + to_string(i));
    seg.getEdge();
    seg.initHoughSpace();
    seg.lineDetect();
    seg.drawLines();
    seg.getCrossPoints();
    seg.warping(img);
  }
  return 0;
}