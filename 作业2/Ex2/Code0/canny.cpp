#include "canny.h"
#include <iostream>
#include <cmath>
#include <stack>
using namespace std;

#define M_PI 3.14159265358979323846
#define BOOSTBLURFACTOR 90.0
#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0

canny::canny(string filename, float _sigma, float _tlow, float _thigh) {
  cout << "Welcome to Canny!" << endl;
  const char* file = filename.c_str();
  read_pgm_image(file);

  sigma = _sigma;
  tlow = _tlow;
  thigh = _thigh;
}
canny::~canny() {
  cout << "Bye!" << endl;
}

int canny::read_pgm_image(const char* inFileName) {
  image.load(inFileName);
  // set image attributes
  rows = image._width;
  columns = image._height;

  // transform to grey image(for bmp image)
  // Gray=R*0.3+G*0.59+B*0.11
  cimg_forXY(image, x, y) {
    int r = (int)image(x, y, 0);
    int g = (int)image(x, y, 1);
    int b = (int)image(x, y, 2);
    double grey = r*0.3 + g * 0.59 + b*0.11;
    image(x, y, 0) = grey;
    image(x, y, 1) = grey;
    image(x, y, 2) = grey;
  }
  
  image.display("original");
}

int canny::write_pgm_image(char* outFileName) {
  image.save(outFileName);
}

double canny::angle_radians(double x, double y) {
  double xu, yu, ang;

   xu = fabs(x);
   yu = fabs(y);

   if ((xu == 0) && (yu == 0)) return 0;

   ang = atan(yu/xu);

   if (x >= 0) {
      if (y >= 0)
        return ang;
      else
        return (2 * M_PI - ang);
   }
   else {
      if (y >= 0)
        return (M_PI - ang);
      else
        return (M_PI + ang);
   }
}

/*
delta_x: The first devivative image, x-direction.
delta_y: The first devivative image, y-direction.
dir_radians: Gradient direction image.
*/
void canny::radian_direction(int xdirtag, int ydirtag) {
  int r, c, pos;
  double dx, dy;

  cimg_forXY(dir_radians, x, y) {
    dx = (double)delta_x(x, y, 0); // The gradient at (x, y), x-direction
    dy = (double)delta_y(x, y, 0); // The gradient at (x, y), y-direction

    if (xdirtag == 1)
      dx = -dx;
    if (ydirtag == -1)
      dy = -dy;

    float gradient = (float)angle_radians(dx, dy);
    dir_radians(x, y, 0) = gradient;
    dir_radians(x, y, 1) = gradient;
    dir_radians(x, y, 2) = gradient;
  }
}

void canny::derrivative_x_y() {
  cout << "Start derivative!" << endl;
  int r, c, pos;
  CImg<short int> temp1(rows, columns);
  CImg<short int> temp2(rows, columns);
  
  // compute x-derivative
  cimg_forXY(temp1, x, y) {
    int x_d;
    if (x == 0) {
      x_d = (int)smoothedim(x+1, y) - (int)smoothedim(x, y);
    }
    else if (x == rows - 1) {
      x_d = (int)smoothedim(x, y) - (int)smoothedim(x-1, y);
    }
    else {
      x_d = (int)smoothedim(x+1, y) - (int)smoothedim(x-1, y);
    }

    temp1(x, y) = x_d;
  }

  // compute y-derivative
  cimg_forXY(temp2, x, y) {
    int y_d;
    if (y == 0) {
      y_d = (int)smoothedim(x, y+1) - (int)smoothedim(x, y);
    }
    else if (y == columns - 1) {
      y_d = (int)smoothedim(x, y) - (int)smoothedim(x, y-1);
    }
    else {
      y_d = (int)smoothedim(x, y+1) - (int)smoothedim(x, y-1);
    }
    
    temp2(x, y) = y_d;
  }

  delta_x = temp1;
  delta_y = temp2;

  cout << "Finish derivative!" << endl;
}

void canny::magnitude_x_y() {
  cout << "Start magnitude_x_y" << endl;
  int r, c, pos, sq1, sq2;
  CImg<short int> temp(rows, columns);

  // compute magnitude at(x, y)
  cimg_forXY(temp, x, y) {
    sq1 = (int)delta_x(x, y) * (int)delta_x(x, y);
    sq2 = (int)delta_y(x, y) * (int)delta_y(x, y);
    short int mag = (short int)(0.5 + sqrt((float)sq1 + (float)sq2));
    temp(x, y) = mag;
  }

  magnitude = temp;
  cout << "End of magnitude_x_y" << endl;
}

void canny::make_gaussian_kernel() {
  cout << "Start kernel" << endl;

  int center, i;
  float x, fx, sum = 0.0;
  windowsize = 1 + 2 * ceil(2.5 * sigma);
  center = windowsize / 2;

  kernel = new float[windowsize];

  for (i = 0; i < windowsize; i++) {
    x = (float)(i - center);
    fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(2 * M_PI));
    kernel[i] = fx;
    sum += fx;
  }

  // Standarlize
  for (i = 0; i < windowsize; i++) {
    kernel[i] /= sum;
  }
  
  for (i = 0; i < windowsize; i++) {
    cout << kernel[i] << " " << endl;
  }
  cout << "Finish Kernel" << endl;
}

void canny::gaussian_smooth() {
  // smoothedim = image.blur(sigma);
  make_gaussian_kernel();
  cout << "Start Smooth" << endl;
  CImg<short int> temp(rows, columns, 1, 1, 0);

  int rr, cc, center;
  float *tempim, dot, sum;
  center = windowsize / 2;
  tempim = new float[rows * columns];
  cimg_forXY(image, x, y) {
    dot = 0.0;
    sum = 0.0;
    for (cc = (-center); cc <= center; cc++) {
      if (((y + cc) >= 0) && ((y + cc) < columns)) {
        dot += (float)image(x, y+cc) * kernel[cc + center];
        sum += kernel[center + cc];
      }
    }
    tempim[x*columns + y] = dot/sum;
  }

  cimg_forXY(image, x, y) {
    sum = 0.0;
    dot = 0.0;
    for(rr = -center; rr <= center; rr++) {
      if (((x + rr) >= 0) && ((x + rr) < rows)) {
        dot += tempim[(x + rr)*columns+y] * kernel[center+rr];
        sum += kernel[center + rr];
      }
    }
    temp(x, y) = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
  }

  delete []tempim;
  smoothedim = temp;
  
  smoothedim.display("smooth");
  smoothedim.save("smooth.pgm");
  cout << "Finish smooth" << endl;
}

void canny::non_max_supp() {
  cout << "Start non_max_supp" <<endl;
  CImg<unsigned char> temp = image;
  // the component of the gradient
  int gx,gy;
  // the temp varialbe
  int g1,g2,g3,g4;
  double weight;
  double dTemp,dTemp1,dTemp2;

  cimg_forXY(temp, x, y) {
    if (x == 0 || y == 0 || x == columns -1 || y == rows - 1) {
      temp(x, y) = 0;
    }
  }

  for (int x = 1; x < rows - 1; x++) {
    for (int y = 1; y < columns - 1; y++) {
      if ((int)magnitude(x, y) == 0) {
        // cout << "mag = 0" << endl;
        temp(x, y, 0) = 255;
        temp(x, y, 1) = 255;
        temp(x, y, 2) = 255;
      }
      else {
        dTemp = magnitude(x, y, 0);
        gx = delta_x(x, y, 0);
        gy = delta_y(x, y, 0);

        if (abs(gy) > abs(gx)) {
          weight = fabs(gx) / fabs(gy);
          g2 = magnitude(x, y-1, 0);
          g4 = magnitude(x, y+1, 0);
          if (gx * gy > 0) {
            g1 = magnitude(x-1,y-1,0);
            g3 = magnitude(x+1,y+1,0);
          }
          else {
            g1 = magnitude(x+1,y-1,0);
            g3 = magnitude(x-1,y+1,0);
          }
        }
        else {
          weight = fabs(gy) / fabs(gx);
          g2 = magnitude(x+1,y,0);
          g4 = magnitude(x-1,y,0);
          if (gx * gy > 0) {
            g1 = magnitude(x+1,y+1,0);
            g3 = magnitude(x-1,y-1,0);
          }
          else {
            g1 = magnitude(x-1,y-1,0);
            g3 = magnitude(x+1,y+1,0);
          }
        }

        dTemp1 = weight * g1 + (1 - weight) * g2;
        dTemp2 = weight * g3 + (1 - weight) * g4;

        if (dTemp >= dTemp1 && dTemp >= dTemp2) {
          temp(x, y, 0) = POSSIBLE_EDGE;
          temp(x, y, 1) = POSSIBLE_EDGE;
          temp(x, y, 2) = POSSIBLE_EDGE;
        }
        else {
          temp(x, y, 0) = NOEDGE;
          temp(x, y, 1) = NOEDGE;
          temp(x, y, 2) = NOEDGE;
        }
      }
    }
  }

  temp.display("nms");
  temp.save("nms.pgm");
  nms = temp;
  cout << "End of non_max_supp" << endl;
}

void canny::follow_edges(CImg<unsigned char> &edgemap, int x, int y, int lowval) {
  int i;
  float thethresh;
  int dx[8] = {1,1,0,-1,-1,-1,0,1},
       dy[8] = {0,1,1,1,0,-1,-1,-1};
  for (i = 0; i < 8; i++) {
    x += dx[i];
    y += dy[i];
    if (((int)edgemap(x, y) == POSSIBLE_EDGE) && ((int)magnitude(x, y) > lowval)) {
      edgemap(x, y) = (unsigned char) EDGE;
      follow_edges(edgemap, x, y, lowval);
    }
  }
}

void canny::apply_hysteresis() {
  cout << "Start Hysteresis" << endl;
  CImg<unsigned char> temp(rows, columns);

  int hist[32768] = {0};
  int r, numedges, maximum_mag, highcount, highthreshold, lowthreshold;
  cimg_forXY(temp, x, y) {
    if ((int)nms(x, y) == POSSIBLE_EDGE) {
      temp(x, y) = POSSIBLE_EDGE;
    }
    else {
      temp(x, y) = NOEDGE;
    }
  }
  edge = temp;

  // handle border
  cimg_forXY(edge, x, y) {
    if (x == 0 || y == 0 || x == rows-1 || y == columns - 1) {
      edge(x, y) = NOEDGE;
    }
  }

  cimg_forXY(edge, x, y) {
    if ((int)edge(x, y) == POSSIBLE_EDGE) {
      hist[magnitude(x, y)]++;
    }
  }

  for (r = 1, numedges = 0; r < 32768; r++) {
    if (hist[r] != 0) {
      maximum_mag = r;
    }
    numedges += hist[r];
  }

  highcount = (int)(numedges * thigh + 0.5);

  r = 1;
  numedges = hist[1];
  while ((r < (maximum_mag-1)) && (numedges < highcount)) {
    r++;
    numedges += hist[r];
  }
  highthreshold = r;
  lowthreshold = (int)(highthreshold * tlow + 0.5);

  cimg_forXY(edge, x, y) {
    if ((int)edge(x, y) == POSSIBLE_EDGE && (int)magnitude(x, y) >= highthreshold) {
      edge(x, y) = EDGE;
      follow_edges(edge, x, y, lowthreshold);
    }
  }

  cimg_forXY(edge, x, y) {
    if (edge(x, y) != EDGE) {
      edge(x, y) = NOEDGE;
    }
  }
  edge.display("hysteresis");
  edge.save("hysteresis.pgm");
  cout << "End of Hysteresis" << endl;
}

void canny::modify() {
  cout << "Start modify" << endl;
  for (int i= 2;i < rows - 2;i++) {
    for (int j = 2;j < columns - 2;j++) {
      //如果该中心点为255,则考虑它的八邻域
      if ((int)edge(i, j) == EDGE) {
        int num = 0;
        for (int k = -1;k < 2;k++) {
          for (int l = -1;l < 2;l++) {
            //如果八邻域中有灰度值为0的点，则去找该点的十六邻域
            if(k != 0 && l != 0 && (int)edge(i+k, j+l) == EDGE)
              num++;
          }
        }

        //如果八邻域中只有一个点是255，说明该中心点为端点，则考虑他的十六邻域
        if(num == 1) {
          for (int k = -2;k < 3;k++) {
            for (int l = -2;l < 3;l++) {
              //如果该点的十六邻域中有255的点，则该点与中心点之间的点置为255
              if (!(k < 2 && k > -2 && l < 2 && l > -2) && (int)edge(i+k, j+l) == EDGE) {
                // cout << "Add edge" << endl;
                edge(i+k/2, j+l/2) = EDGE;
              }
            }
          }
        }
      }
    }
  }
  cout << "Link!" << endl;
  edge.display("link");


  /*eliminate edges smaller than 20*/
  /*
  vector<vector<bool>> isvisited(rows, vector<bool>(columns));
  vector<stack<Point>> path(1500, stack<Point>());
  queue<Point> q;
  int count = 0;

  // Initialize the flag matrix
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      if ((int)edge(i, j) == NOEDGE) {
        isvisited[i][j] = true;
      }
      else {
        isvisited[i][j] = false;
      }
    }
  }

  int dx[8] = {1,1,0,-1,-1,-1,0,1},
       dy[8] = {0,1,1,1,0,-1,-1,-1};

  cout << "Start bfs" << endl;
  // start bfs
  int pos = 0;
  cimg_forXY(edge, x, y) {
    count = 0;
    if (isvisited[x][y] == false) {
      q.push(Point(x, y));
      path[pos].push(Point(x, y));
      while (!q.empty()) {
        Point u = q.front();
        q.pop();
        if (!isvisited[u.x][u.y]) {
          isvisited[u.x][u.y] = true;
          path[pos].push(u);
          count++;
          for (int i = 0; i < 8; i++) {
            int x_d = u.x + dx[i];
            int y_d = u.y + dy[i];
            if (x_d >= 0 && x_d < rows && y_d >= 0 && y_d < columns) {
              if ((isvisited[x_d][y_d] == false)) {
                q.push(Point(x_d, y_d));
              }
            }
          }
        }
      }
      if (path[pos].size() <= 500) {
        while (!path[pos].empty()) {
          Point temp = path[pos].top();
          edge(temp.x, temp.y) = NOEDGE;
          path[pos].pop();
        }
      }
      pos++;
    }
  }
  edge.display("final");
  */
  edge.save("final.bmp");
  cout << "End of modify" << endl;
}


int main(int argc,char *argv[]) {
  if (argc < 5) {
    cout << "\n<USAGE> " << argv[0] <<  " image sigma tlow thigh [writedirim]\n";
    cout << "\n      image:      An image to process. Must be in BMP format.\n";
    cout << "      sigma:      Standard deviation of the gaussian blur kernel.\n";
    cout << "      tlow:       Fraction (0.0-1.0) of the high edge strength threshold.\n";
    cout << "      thigh:      Fraction (0.0-1.0) of the distribution of non-zero edge\n";
    cout << "                  strengths for hysteresis. The fraction is used to compute\n";
    cout << "                  the high edge strength threshold.\n";
    return 0;
  }

  string infilename = argv[1];
  float sigma = atof(argv[2]);
  float tlow = atof(argv[3]);
  float thigh = atof(argv[4]);

  canny c1(infilename, sigma, tlow, thigh);
  c1.gaussian_smooth();
  c1.derrivative_x_y();
  c1.radian_direction(-1, -1);
  c1.magnitude_x_y();
  c1.non_max_supp();
  c1.apply_hysteresis();
  c1.modify();
  return 0;
}