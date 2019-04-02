#include <iostream>
#include "hw1.hpp"
using namespace std;

int main() {
  // Prompt user to choose
  int cmd;
  cout << "Please input your choice:" << endl;
  cout << "[0]: use CImg functions" << endl;
  cout << "[1]: without using CImg functions" << endl;
  cout << "Your Choice: " << endl;
  cin >> cmd;

  MyImage m;
  if (cmd == 0) {
    m.changeColor();
    m.drawBlueCircle();
    m.drawYellowCircle();
    m.drawBlueLine();
  }
  else if (cmd == 1) {
    m.changeColor();
    m.drawBlueCircleByMe();
    m.drawYellowCircleByMe();
    m.Bresenham();
  }
  cout << "Please check 2.bmp in your file, Thank you!" << endl;
  return 0;
}