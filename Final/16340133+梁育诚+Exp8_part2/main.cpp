#include <iostream>
#include <vector>
#include <list>
#include <iostream>
#include <stdio.h>
#include "Segmentation.h"
#include "ImageSplit.h"
using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    return 0;
  }
  int i = atoi(argv[1]);
  string root = "ImageData/";
  string warpingPath = "./result/warping/";
  //for (int i = 10; i <= 10; i++) {
    CImg<float> img((root + to_string(i) + ".bmp").c_str());
    CImg<float> grayImg(img._width, img._height, 1, 1);
    
    // Change to Grey Image
    cimg_forXY(img, x, y) {
      grayImg(x, y, 0, 0) = 0.299 * img(x, y, 0, 0) + 0.587 * img(x, y, 0, 1) + 0.114 * img(x, y, 0, 2);
    }

    // 矫正A4纸
    Segmentation seg(grayImg);
    seg.segment();
    seg.binarization("Seg" + to_string(i));
    seg.getEdge();
    seg.initHoughSpace();
    seg.lineDetect();
    seg.drawLines();
    seg.getCrossPoints();
    CImg<float> a4 = seg.warping(img);
    //a4.display();
    a4.save((warpingPath + to_string(i) + ".bmp").c_str());
  
    // 提取A4纸中的数字
    string singleNumPath = "result/singleNumImg/";
    singleNumPath = singleNumPath + to_string(i) + ".bmp";
    CImg<float> image((warpingPath + to_string(i) + ".bmp").c_str());
    ImageSplit splitImage;
    splitImage.run(a4, singleNumPath.c_str());
  //}
  return 0;
}