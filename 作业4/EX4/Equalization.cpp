#include "Equalization.hpp"
#define E 2.71828
#define PI 3.14159265358979323846

Equalization::Equalization(string so) {
  // Read image
  const char * s = so.c_str();
  image.load_bmp(s);
  imageRows = image._width;
  imageColumns = image._height;
  // Display
  image.display("image");

  // To grey
  
  cimg_forXY(image, x, y) {
    int R = image(x, y, 0);
    int G = image(x, y, 1);
    int B = image(x, y, 2);

    double grey = R * 0.2126 + G * 0.7152 + B * 0.0722;
    image(x, y, 0) = grey;
    image(x, y, 1) = grey;
    image(x, y, 2) = grey;
  }
  // Display grey image
  image.display("grey");
  
}

Equalization::Equalization(string so, string tar) {
  // Read source
  const char *s = so.c_str();
  origin.load_bmp(s);
  originRows = origin._width;
  originColumns = origin._height;

  // Read target
  const char *t = tar.c_str();
  target.load_bmp(t);
  targetRows = target._width;
  targetColumns = target._height;

  result = origin;

  // Display
  origin.display("origin", false);
  target.display("target", false);
}

Equalization::~Equalization() {

}


void Equalization::equalization() {
  float statistics_R[256] = {0};
  float cumulative_R[256] = {0};
  float probability_R[256] = {0};

  float statistics_G[256] = {0};
  float cumulative_G[256] = {0};
  float probability_G[256] = {0};

  float statistics_B[256] = {0};
  float cumulative_B[256] = {0};
  float probability_B[256] = {0};

  // Display Origin histogram
  CImg<float> temp = image;
  CImg<int> hist = temp.histogram(256, 0, 255);
  hist.display_graph("Origin Histogram", 3);
  
  output = image;

  // Total Pixel
  int total = imageRows * imageColumns;

  // Compute Statistic
  cimg_forXY(image, x, y) {
    int R = (int)image(x, y, 0);
    int G = (int)image(x, y, 1);
    int B = (int)image(x, y, 2);
    //cout << R << " " << G << " " << B << endl;
    statistics_R[R]++;
    statistics_G[G]++;
    statistics_B[B]++;
  }

  // Compute Probability
  for (int i = 0; i < 256; i++) {
    probability_R[i] = statistics_R[i] / total;
    probability_G[i] = statistics_G[i] / total;
    probability_B[i] = statistics_B[i] / total;
  }

  // Compute cumulative
  for (int i = 0; i < 256; i++) {
    for (int j = 0; j < i; j++) {
      cumulative_R[i] += probability_R[j];
      cumulative_G[i] += probability_G[j];
      cumulative_B[i] += probability_B[j];
    }
  }

  // Color Transformation
  cimg_forXY(output, x, y) {
    int R = image(x, y, 0);
    int G = image(x, y, 1);
    int B = image(x, y, 2);
    output(x, y, 0) = 255 * cumulative_R[R];
    output(x, y, 1) = 255 * cumulative_G[G];
    output(x, y, 2) = 255 * cumulative_B[B];
  }

  // Display Output histogram
  temp = output;
  hist = temp.histogram(256, 0, 255);
  hist.display_graph("Output Histogram", 3);

  output.display("output", false);
}

void Equalization::colorTransform() {
  float sum1[3], sum2[3], sum1_squ[3], sum2_squ[3];
  float mean1[3], mean2[3],vari1[3], vari2[3];

  // Total Pixels
  int total1 = originRows * originColumns;
  int total2 = targetRows * targetColumns;

  // Initialize variables
  for (int i = 0; i < 3; i++) {
    sum1[i] = sum2[i] = sum1_squ[i] = sum2_squ[i] = 0;
    mean1[i] = mean2[i] = vari1[i] = vari2[i] = 0;
  }

  cout << "Origin: RGB to LAB" << endl;
  // RGB to LAB(origin)
  cimg_forXY(origin, x, y) {
    float R = origin(x, y, 0);
    float G = origin(x, y, 1);
    float B = origin(x, y, 2);

    float L = 0.3811 * R + 0.5783 * G + 0.0402 * B;
    float M = 0.1967 * R + 0.7244 * G + 0.0782 * B;
    float S = 0.0241 * R + 0.1288 * G + 0.8444 * B;

    if (L == 0) {
      L = 1;
    }
    if (M == 0) {
      M = 1;
    }
    if (S == 0) {
      S = 1;
    }

    L = log(L);
    M = log(M);
    S = log(S);

    float l = (1.0 / sqrt(3)) * (L + M + S);
    float a = 1.0 / sqrt(6) * L + 1.0 / sqrt(6) * M - 2.0 / sqrt(6) * S;
    float b = 1.0 / sqrt(2) * L - 1.0 / sqrt(2) * M;

    result(x, y, 0) = l;
    result(x, y, 1) = a;
    result(x, y, 2) = b;

    // Calculation total values for each channel
    sum1[0] += l;
    sum1[1] += a;
    sum1[2] += b;
  }

  // Calculate mean for each channel in origin
  mean1[0] = sum1[0] / total1;
  mean1[1] = sum1[1] / total1;
  mean1[2] = sum1[2] / total1;

  // Calculate Variance
  cimg_forXY(result, x, y) {
    sum1_squ[0] += (result(x, y, 0) - mean1[0]) * (result(x, y, 0) - mean1[0]);
    sum1_squ[1] += (result(x, y, 1) - mean1[1]) * (result(x, y, 1) - mean1[1]);
    sum1_squ[2] += (result(x, y, 2) - mean1[2]) * (result(x, y, 2) - mean1[2]);
  }

  vari1[0] = sqrt(sum1_squ[0] / total1);
  vari1[1] = sqrt(sum1_squ[1] / total1);
  vari1[2] = sqrt(sum1_squ[2] / total1);


  cout << "target: RGB to LAB" << endl;
  // RGB to LAB(target)
  cimg_forXY(target, x, y) {
    float R = target(x, y, 0);
    float G = target(x, y, 1);
    float B = target(x, y, 2);

    float L = 0.3811 * R + 0.5783 * G + 0.0402 * B;
    float M = 0.1967 * R + 0.7244 * G + 0.0782 * B;
    float S = 0.0241 * R + 0.1288 * G + 0.8444 * B;

    if (L == 0) {
      L = 1;
    }
    if (M == 0) {
      M = 1;
    }
    if (S == 0) {
      S = 1;
    }

    L = log(L);
    M = log(M);
    S = log(S);

    float l = (1.0 / sqrt(3)) * (L + M + S);
    float a = 1.0 / sqrt(6) * L + 1.0 / sqrt(6) * M - 2.0 / sqrt(6) * S;
    float b = 1.0 / sqrt(2) * L - 1.0 / sqrt(2) * M;

    target(x, y, 0) = l;
    target(x, y, 1) = a;
    target(x, y, 2) = b;

    // Calculation total values for each channel
    sum2[0] += l;
    sum2[1] += a;
    sum2[2] += b;
  }

  // Calculate mean for each channel in target
  mean2[0] = sum2[0] / total2;
  mean2[1] = sum2[1] / total2;
  mean2[2] = sum2[2] / total2;

  // Calculate Variance
  cimg_forXY(target, x, y) {
    sum2_squ[0] += (target(x, y, 0) - mean2[0]) * (target(x, y, 0) - mean2[0]);
    sum2_squ[1] += (target(x, y, 1) - mean2[1]) * (target(x, y, 1) - mean2[1]);
    sum2_squ[2] += (target(x, y, 2) - mean2[2]) * (target(x, y, 2) - mean2[2]);
  }

  vari2[0] = sqrt(sum2_squ[0] / total2);
  vari2[1] = sqrt(sum2_squ[1] / total2);
  vari2[2] = sqrt(sum2_squ[2] / total2);


  cout << "Color Transformation"<< endl;
  // Color Transformation!
  cimg_forXY(result, x, y) {
    result(x, y, 0) = (result(x, y, 0) - mean1[0]) * vari2[0] / vari1[0] + mean2[0];
    result(x, y, 1) = (result(x, y, 1) - mean1[1]) * vari2[1] / vari1[1] + mean2[1];
    result(x, y, 2) = (result(x, y, 2) - mean1[2]) * vari2[2] / vari1[2] + mean2[2];
  }

  // LAB to RGB(result)
  cimg_forXY(result, x, y) {
    float L = (1.0 / sqrt(3)) * result(x, y, 0) + (1.0 / sqrt(6)) * result(x, y, 1) 
              + (1.0 / sqrt(2)) * result(x, y, 2);
    float M = (1.0 / sqrt(3)) * result(x, y, 0) + (1.0 / sqrt(6)) * result(x, y, 1) 
              - (1.0 / sqrt(2)) * result(x, y, 2);
    float S = (1.0 / sqrt(3)) * result(x, y, 0) - (sqrt(6) / 3) * result(x, y, 1);

    L = pow(E, L);
    M = pow(E, M);
    S = pow(E, S);

    result(x, y, 0) = 4.4679 * L - 3.5873 * M + 0.1193 * S;
    result(x, y, 1) = -1.2186 * L + 2.3809 * M - 0.1624 * S;
    result(x, y, 2) = 0.0497 * L - 0.2439 * M + 1.2045 * S;
  }

  result.display("result", false);
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    return 0;
  }
  string source = argv[1];
  string target = argv[2];
  int op = atoi(argv[3]);
  if (op == 0) {
    Equalization obj(source);
    obj.equalization();
  }
  else {
    Equalization obj(source, target);
  obj.colorTransform();
  }
  
  return 0;
}