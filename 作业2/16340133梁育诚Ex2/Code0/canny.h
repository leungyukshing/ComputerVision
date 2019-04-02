#include <string>
#include "CImg.h"
#include <queue>

using namespace cimg_library;
using namespace std;

struct Point {
  int x;
  int y;
  Point(int _x, int _y) {
    x = _x;
    y = _y;
  }
};

class canny
{
public:
  canny(string, float, float, float);
  ~canny();

  int read_pgm_image(const char* inFileName);
  int write_pgm_image(char* outFieName);
  double angle_radians(double x, double y);
  void radian_direction(int xdirtag, int ydirtag);
  void derrivative_x_y();
  void magnitude_x_y();
  void make_gaussian_kernel();
  void gaussian_smooth();
  void non_max_supp();
  void follow_edges(CImg<unsigned char> &edgemap, int x, int y, int lowval);
  void apply_hysteresis();
  void modify();
private:
  string inFileName;
  char* dirFileName;
  char* outFileName;
  char* composedFName;
  CImg<unsigned char> image;
  CImg<unsigned char> edge;
  CImg<float> dirim; // store the direction of gradient
  int rows;
  int columns;
  float sigma,
    tlow,
    thigh;

  CImg<short int> delta_x; /* The first devivative image, x-direction. */
  CImg<short int> delta_y; /* The first derivative image, y-direction. */
  CImg<float> dir_radians; /* Gradient direction image.                */
  CImg<short int> smoothedim; /* The image after gaussian smoothing.      */
  CImg<short int> magnitude; /* The magnitude of the gadient image.      */
  float *kernel; /* A two dimensional gaussian kernel. */
  CImg<unsigned char> nms; /* Points that are local maximal magnitude. */
  int windowsize;/* Dimension of the gaussian kernel. */
};