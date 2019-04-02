#include "Segmentation.h"

Segmentation::Segmentation(CImg<float> &_img) {
  img = _img;
  rows = img._width;
  columns = img._height;

  x_max = img._width - 1;
  y_max = img._height - 1;
  gradient_threshold = 0;
  vote_threshold = 300;
  peak_dis = 20;

  // Initialize histogram
  for (int i = 0; i < 256; i++) {
    histogram[i] = 0;
  }

  // Initialize sin and cos sets
  for (int i = 0; i < 360; i++) {
    sins.push_back(sin(2 * M_PI * i / 360));
    coss.push_back(cos(2 * M_PI * i / 360));
  }
}

float Segmentation::getThreshold() { return threshold; }

// Compute the histogram of the image
void Segmentation::getHistogram() {
  cimg_forXY(img, x, y) {
    histogram[(int)img(x, y)]++;
  }
}

// Do segmentation
void Segmentation::segment() {
  int totalPixels = img._width * img._height;
  int gray = 0, grayBackground = 0, grayFrontground = 0;
  float wBackground = 0, wFrontground = 0;
  float g_best = 0;
  float g = 0;
  int numBackground = 0, numFrontground = 0;
    
  getHistogram();
  for (int i = 0; i < 256; i++) {
    gray += i * histogram[i];
  }

  for (int i = 0; i < 256; i++) {
    grayBackground += histogram[i] * i; // 背景总灰度值
    grayFrontground = gray - grayBackground; // 前景总灰度值
    numBackground += histogram[i]; // 背景pixel数量
    numFrontground = totalPixels - numBackground; // 前景pixel数量
    wBackground = numBackground * 1.0 / totalPixels; // 背景pixel占比
    wFrontground = 1 - wBackground; // 前景pixel占比
      
    // 无法区分前景和背景
    if (wBackground == 0 || wFrontground == 0)
      continue;
    // 找出最佳阈值
    g = wBackground * wFrontground * (grayBackground / numBackground - grayFrontground / numFrontground) 
          * (grayBackground / numBackground - grayFrontground / numFrontground);
    if (g > g_best) {
      threshold = i;
      g_best = g;
    }
  }
  cout << "Threshold: " << threshold << endl;
}

// Detect Edge
void Segmentation::getEdge() {
  
  // using sobel to compute gradient
  CImg<float> sobelx(3, 3, 1, 1, 0);
  CImg<float> sobely(3, 3, 1, 1, 0);
  // Initialize sobel
  sobelx(0, 0) = -1, sobely(0, 0) = 1;
  sobelx(0, 1) = 0, sobely(0, 1) = 2;
  sobelx(0, 2) = 1, sobely(0, 2) = 1;
  sobelx(1, 0) = -2, sobely(1, 0) = 0;
  sobelx(1, 1) = 0, sobely(1, 1) = 0;
  sobelx(1, 2) = 2, sobely(1, 2) = 0;
  sobelx(2, 0) = -1, sobely(2, 0) = -1;
  sobelx(2, 1) = 0, sobely(2, 1) = -2;
  sobelx(2, 2) = 1, sobely(2, 2) = -1;

  CImg<float> gradient_x = result;
  gradient_x = gradient_x.get_convolve(sobelx);
  CImg<float> gradient_y = result;
  gradient_y = gradient_y.get_convolve(sobely);

  cimg_forXY(result, i, j) {
    double grad = sqrt(gradient_x(i, j) * gradient_x(i, j) + gradient_y(i, j) * gradient_y(i, j));
    //cout << "Grad: " << grad << endl;
    if (grad > gradient_threshold) {
      result(i, j) = EDGE;
    }
    else {
      result(i, j) = NOEDGE;
    }
  }
  
  /*
  int rows = result._width;
  int cols = result._height;
  vector<vector<int>> detX(rows, vector<int>(cols, 0));

    /* Compute the x-derivative filter = [-1, 0, +1] 
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (c == 0) detX[r][c] = result(r, c+1) - result(r, c);
            else if (c == cols-1) detX[r][c] = result(r, c) - result(r, c-1);
            else detX[r][c] = result(r, c+1) - result(r, c-1);
        }
    }

    vector<vector<int>> detY(rows, vector<int>(cols, 0));

    /* Compute the y-derivative filter = [-1, 0, +1]^T 
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (r == 0) detY[r][c] = result(r+1, c) - result(r, c);
            else if (r == rows-1) detY[r][c] = result(r, c) - result(r-1, c);
            else detY[r][c] = result(r+1, c) - result(r-1, c);
        }
    }    

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (detX[r][c] != 0 || detY[r][c] != 0) 
                result(r, c) = 0;
            else
                result(r, c) = 255;
        }
    }
    */
  //result.display("edge", false);
}

// Initialize Hough Space
void Segmentation::initHoughSpace() {
  int maxp = (int)sqrt(pow(result._width / 2, 2) + pow(result._height / 2, 2));
  hough_space = CImg<float>(maxp, 360);
  hough_space.fill(0);

  cimg_forXY(result, x, y) {
    int value = (int)result(x, y);
    int p; 
    // Only care the edge points
    if (value == EDGE) {
      // transform to polar coordinate (p, i) => (length，angle)
      int x0 = x - rows / 2;
      int y0 = columns / 2 - y;
      // Calculate all lines gone through a point(all angles)
      for (int i = 0; i < 360; i++) {
        p = x0 * coss[i] + y0 * sins[i];
        // Voting
        if (p >= 0 && p < maxp) {
          hough_space(p, i)++;
        }
      }
    }
  }
}

// Detect Line from houghspace
void Segmentation::lineDetect() {
  cout << "Begin detect" << endl;
  int wid = hough_space._width;
  int hei = hough_space._height;
  int maxVotes;
  for (int i = 0; i < wid; i += PART / 2) {
    for (int j = 0; j < hei; j += PART / 2) {
      maxVotes = getPartMax(hough_space, i, j);
      for (int x = i; (x < i + PART) && (x < wid); x++) {
        for (int y = j; (y < j + PART) && (y < hei); y++) {
          if ((int)hough_space(x, y) < maxVotes) {
            hough_space(x, y) = 0;
          }
        }
      }
    }
  }
  cout << "finish detect" << endl;
  //hough_space.display("detected", false);
  // save the bright points in houghspace
  cimg_forXY(hough_space, x, y) {
    if ((int)hough_space(x, y) != 0) {
      lines.push_back(make_pair(x, y));
      lineCount.push_back((int)hough_space(x, y));
    }
  }
}

// Find local maximum
int Segmentation::getPartMax(CImg<float>& img, int &x, int &y) {
  int wid = (x + PART > img._width) ? img._width : x + PART;
  int hei = (y + PART > img._height) ? img._height : y + PART;
  int maxVotes = 0;
  for (int i = x; i < wid; i++) {
    for (int j = y; j < hei; j++) {
      maxVotes = ((int)img(i, j) > maxVotes) ? ((int)img(i, j)) : maxVotes;
    }
  }
  return maxVotes;
}

// draw the outline of A4
void Segmentation::drawLines() {
  cout << "Begin drawLines" << endl;
  int maxLength;
  maxLength = sqrt(pow(rows / 2, 2) + pow(columns / 2, 2)); // width of houghspace

  showEdge = CImg<float>(rows, columns, 1, 1, 0); // Initialize showEdge
  afterHough = CImg<float>(rows, columns, 1, 3, 0); // Initialize output image with original image
    
  sortedLineCount = lineCount; // Initialize the vector

  // Sort lines by its votes, descending
  sort(sortedLineCount.begin(), sortedLineCount.end(), greater<int>());

  // Show the number of decteced size
  cout << "Size: " << sortedLineCount.size() << endl;

  // Record the parameters of the lines
  vector<pair<int, int>> points;
  for (int i = 0; i < 4; i++) {
    int weight = sortedLineCount[i];
    int indexInLines; // Index in lines vector
    vector<int>::iterator it;
    it = find(lineCount.begin(), lineCount.end(), weight);
    indexInLines = it - lineCount.begin();
    points.push_back(lines[indexInLines]);
  }
    
  // Draw lines in output image
  cout << endl << endl;
  for (int i = 0; i < points.size(); i++) {
    int p = points[i].first;
    int angle = points[i].second;
    double k;
    double b;
    if (sins[angle] == 0) {
      b = p / coss[angle];
      cout << "Line " << i << ": x = " << b << endl;
    }
    else {
      k = -(coss[angle] * 1.0 / sins[angle]);
      b = p * 1.0 / sins[angle];
      cout << "Line " << i << ": y = " << k << " * x + " << b << endl;
    }
    
    cout << afterHough._width << " " << afterHough._height << endl;
    // dye the points that satisfy the equation with BLUE
    cimg_forXY(showEdge, x, y) {
      int x0 = x - rows / 2;
      int y0 = columns / 2 - y;
      int p_ = x0 * coss[angle] + y0 * sins[angle];
      
      if (p == p_) {
        showEdge(x, y) += 255.0 / 2;
        // line is BLUE
        afterHough(x, y, 0) = 0;
        afterHough(x, y, 1) = 0;
        afterHough(x, y, 2) = 255;
      }
    }
  }
  cout << endl << endl;
  //afterHough.display("afterHough", false);
  cout << "End of drawLines" << endl;
}

// Get crosspoint
void Segmentation::getCrossPoints() {
  cout << "Begin drawPoints" << endl;
  unsigned char green[3] = {0, 255, 0};
  for (int x = 0; x < rows-1; x++) {
    for (int y = 0; y < columns-1; y++) {
      int area[4];
      area[0] = (int)showEdge(x, y);
      area[1] = (int)showEdge(x + 1, y);
      area[2] = (int)showEdge(x, y + 1);
      area[3] = (int)showEdge(x + 1, y + 1);
      // Enough point, indicates a crosspoint!
      if (area[0] + area[1] + area[2] + area[3] >= 255*3/2) {
        bool flag = true;
        for (int i = 0; i < intersections.size(); i++) {
          double dis = pow((x-intersections[i].x), 2) + pow((y-intersections[i].y), 2);
          dis = sqrt(dis);
          if (dis <= 10) {
            flag = false;
            break;
          }
        }
        if (flag) {
          intersections.push_back(Point(x, y, 0));
        }
        afterHough.draw_circle(x, y, 5, green);
      }
    }
  }
  cout << "End of drawPoints" << endl;
  //afterHough.display("drawPoints", false);

  // 对四个角点排序
  // 先对y排序
  sort(intersections.begin(), intersections.end());
  // 对x分别排序
  if (intersections[0].x > intersections[1].x) {
    Point p = intersections[0];
    intersections[0] = intersections[1];
    intersections[1] = p;
  }
  if (intersections[2].x > intersections[3].x) {
    Point p = intersections[2];
    intersections[2] = intersections[3];
    intersections[3] = p;
  }
  cout << "Points: " << endl;
  for (int i = 0; i < intersections.size(); i++) {
    cout << intersections[i].x << " " << intersections[i].y << endl;
  }
  
}

// 计算透视变换矩阵
vector<CImg<float>> Segmentation::computeTransformMatrix(CImg<float> a4) {
  vector<Point> destRectPoints;
  Point tempPoint1(0, 0, 0);
  Point tempPoint2(a4._width - 1, 0, 0);
  Point tempPoint3(0, a4._height - 1, 0);
  Point tempPoint4(a4._width - 1, a4._height - 1, 0);
  destRectPoints.push_back(tempPoint1);
  destRectPoints.push_back(tempPoint2);
  destRectPoints.push_back(tempPoint3);
  destRectPoints.push_back(tempPoint4);

  CImg<float> y1(1, 3, 1, 1, 0), y2(1, 3, 1, 1, 0), y3(1, 3, 1, 1, 0), y4(1, 3, 1, 1, 0);
  CImg<float> c1(1, 3, 1, 1, 0), c2(1, 3, 1, 1, 0), c3(1, 3, 1, 1, 0), c4(1, 3, 1, 1, 0);
  CImg<float> A1(3, 3, 1, 1, 1), A2(3, 3, 1, 1, 1);

  for (int i = 0; i < 3; i++) {
    A1(0, i) = destRectPoints[i].x; A1(1, i) = destRectPoints[i].y;
    A2(0, i) = destRectPoints[3-i].x; A2(1, i) = destRectPoints[3-i].y;

    y1(0, i) = intersections[i].x; y2(0, i) = intersections[i].y;
    y3(0, i) = intersections[3-i].x; y4(0, i) = intersections[3-i].y;
  }
  c1 = y1.solve(A1); c2 = y2.solve(A1);
  c3 = y3.solve(A2); c4 = y4.solve(A2);

  CImg<float> temptransform1(3, 3, 1, 1, 0), temptransform2(3, 3, 1, 1, 0);
  for (int i = 0; i < 3; i++) {
    temptransform1(i, 0) = c1(0, i);
    temptransform1(i, 1) = c2(0, i);

    temptransform2(i, 0) = c3(0, i);
    temptransform2(i, 1) = c4(0, i);
  }
  temptransform1(0, 2) = 0; temptransform1(1, 2) = 0; temptransform1(2, 2) = 1;
  temptransform2(0, 2) = 0; temptransform2(1, 2) = 0; temptransform2(2, 2) = 1;
  vector<CImg<float> > temptransform;
  temptransform.push_back(temptransform1);
  temptransform.push_back(temptransform2);
  return temptransform;
}

// Perspective Transformation
CImg<float> Segmentation::warping(CImg<float> srcImg) {
  a4 = CImg<float>(srcImg._width, srcImg._height, 1, 3, 0);
    
  vector<CImg<float> > transform;
  transform = computeTransformMatrix(a4); // 计算变换矩阵

  CImg<float> y(1, 2, 1, 1, 0);
  CImg<float> c(1, 2, 1, 1, 0);
  CImg<float> A(2, 2, 1, 1, 1);
  A(0, 0) = 0;
  A(0, 1) = a4._width - 1;
  y(0, 0) = a4._height - 1;
  y(0, 1) = 0;
  c = y.solve(A);

  CImg<float> temp1(1, 3, 1, 1, 1), temp2(1, 3, 1, 1, 1);
  cimg_forXY(a4, i, j) {
    temp1(0, 0) = i;
    temp1(0, 1) = j;
      
    double inner_procuct = i * c(0, 0) - j + c(0, 1);
    temp2 = inner_procuct >= 0 ? transform[0] * temp1 : transform[1] * temp1;
    temp2(0, 0) = temp2(0, 0) < 0 ? 0 : (temp2(0, 0) > x_max ? x_max : temp2(0, 0));
    temp2(0, 1) = temp2(0, 1) < 0 ? 0 : (temp2(0, 1) > y_max ? y_max : temp2(0, 1));

    a4(i, j, 0) = srcImg(temp2(0, 0), temp2(0, 1), 0);
    a4(i, j, 1) = srcImg(temp2(0, 0), temp2(0, 1), 1);
    a4(i, j, 2) = srcImg(temp2(0, 0), temp2(0, 1), 2);
  }
  //a4.display("a4", false);
  return a4;
}

// Get Binary Image
void Segmentation::binarization(string path) {
  result = CImg<float>(img._width, img._height, 1);
  threshold = threshold * 1.2;
  cimg_forXY(img, x, y) {
    if (img(x, y) > threshold) {
      result(x, y, 0, 0) = 0;
    }
    else {
      result(x, y, 0, 0) = 255;
    }
  }
  //result.display("seg");
}