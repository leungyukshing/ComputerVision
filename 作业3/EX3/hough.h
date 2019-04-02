#include <string>
#include "CImg.h"
#include <queue>

using namespace cimg_library;
using namespace std;

struct Point {
  int x;
  int y;
  int z;
  int mag;
  Point(int _x, int _y, int _z, int _mag) {
    x = _x;
    y = _y;
    z = _z;
    mag = _mag;
  }
  bool operator<(const Point & r) const {
    return mag < r.mag;
  }
};

class canny
{
public:
  canny(string, float, float, float, int, int);
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


  // Hough Line Transformation
  void houghLine(); // 转换为霍夫空间
  void lineDetect(); // 在霍夫空间中检测出直线
  int getPartMax(CImg<float>& img, int &x, int &y); // 局部极大值点（直线交点）
  void drawLines();
  void drawPoints();

  // Hough Circle Transformation
  void houghCircle();
  void circleDetect();
  //void drawCircles();
  //void drawCircles(int radius);
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

  // Hough
  vector<double> sins; // value of all sin
  vector<double> coss; // value of all cos
  CImg<float> houghspace; // Hough space image
  CImg<unsigned char> afterHough; // after using hough transformation
  vector<pair<int, int>> lines; // set of points
  vector<int> lineCount; // lines' votes set
  vector<int> sortedLineCount; // sorted lines set by its votes
  CImg<float> showEdge; // image to show the lines detected
  CImg<float> output;

  vector<Point> circles; // Set of circle centers
  int circleNumber;
  int minRadius, maxRadius; // thresholds to filters some circles
};