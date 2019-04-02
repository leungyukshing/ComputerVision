#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "CImg.h"
#include "AlignFeature.h"
#include "Blend.h"
#include "Stitch.h"
#include "ImageWrap.h"


using namespace std;
using namespace cimg_library;

int main() {
  cout << "Image Stitching!" << endl;
  // Input Image Resource
  int setNumber = 0;
  cout << "Please Input Dataset: ";
  cin >> setNumber;

  // DataSet 1
  if (setNumber == 1) {
    string path = "./TEST-ImageData(1)/pano1_000";
    vector<CImg<float> > imgs;
    for (int i = 1; i <= 4; i++) {
      CImg<float> img((path + to_string(i) + ".bmp").c_str());
      imgs.push_back(img);
    }
    CImg<float> result = stitching(imgs);
    result.display("result");

    result.save("./result/result1.bmp");
  }
  // DataSet 2
  else if (setNumber == 2) {
    string path = "./TEST-ImageData(2)/";
    vector<CImg<float> > imgs;
    for (int i = 8; i <= 25; i++) {
      CImg<float> img((path + to_string(i) + ".bmp").c_str());
      imgs.push_back(img);
    }
    CImg<float> result = stitching(imgs);
    result.display("result");
    result.save("./result/result2.bmp");
  }
  // DataSet 3
  /*
  else if (setNumber == 3) {
    string path = "./TEST-ImageData(3)/";
    vector<CImg<float> > imgs;
    for (int i = 1; i <= 5; i++) {
      CImg<float> img((path + to_string(i) + ".bmp").c_str());
      imgs.push_back(img);
    }
    CImg<float> result = stitching(imgs);
    result.display("result");
    result.save("./result/res3.bmp");
  }
  */
  else {
    cout << "DataSet Not Existed!" << endl;
  }
  return 0;
}