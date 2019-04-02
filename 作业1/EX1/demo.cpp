#include "CImg.h"
using namespace cimg_library;
using namespace std;

void test_cimg() {
  // read BMP image
  CImg<unsigned char> SrcImg;
  SrcImg.load_bmp("1.bmp");

  // get width and height
  int w = SrcImg._width;
  int h = SrcImg._height;

  // Initialize a new grey image
  CImg<unsigned char> TempImg(w, h, 1, 1, 0);

  // show image
  SrcImg.display();

  // handle every pixel
  cimg_forXY(SrcImg, x, y) {
    if (SrcImg(x, y, 0) == 102 && SrcImg(x, y, 1) == 102 && SrcImg(x, y, 2) == 102)
      TempImg(x, y) = 255;
  }

  // show new image
  TempImg.display();

  // save image
  TempImg.save("Thred.bmp");
}

int main() {
  test_cimg();
}