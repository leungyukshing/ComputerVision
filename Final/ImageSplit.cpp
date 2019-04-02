#include "ImageSplit.h"
#include <direct.h>

int lineColor[3]{ 255, 0, 0 };

void ImageSplit::run(CImg<float> warpingResult, const string baseAddress) {
  /*
  if (_access(baseAddress.c_str(), 0) == -1)
    _mkdir(baseAddress.c_str());
  */
  basePath = baseAddress + "/"; // 数字分割图片存储地址
  // 提取边缘
  edgeImg = getEdge(warpingResult);
  //edgeImg.display("Edge");

  // 二值化处理
  binaryImg = convertToBinaryImg(edgeImg);
  //binaryImg.display("Binary");

  // 行分割，按照行划分数字
  findDividingLine();
  divideIntoBarItemImg();
  // histogramImg.display("Histogram");
  // dividingImg.display("Divide");

  // 对每张子图操作
  for (int i = 0; i < subImageSet.size(); i++) {
    dilateImg(i); // 对分割后每一张数字做扩张
    connectedRegionsTagging(i); // 连通区域标记算法
    saveSingleNumImg(i); // 存储结果
  }
}

// Detect Edge
CImg<float> ImageSplit::getEdge(CImg<float> img) {
  edgeImg = CImg<float>(img._width, img._height, 1, 1, 0);
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

  CImg<float> gradient_x = img;
  gradient_x = gradient_x.get_convolve(sobelx);
  CImg<float> gradient_y = img;
  gradient_y = gradient_y.get_convolve(sobely);

  cimg_forXY(img, i, j) {
    double grad = sqrt(gradient_x(i, j) * gradient_x(i, j) + gradient_y(i, j) * gradient_y(i, j));
    if (grad > 110) {
      img(i, j, 0) = 0;
      img(i, j, 1) = 0;
      img(i, j, 2) = 0;
    }
    else {
      img(i, j, 0) = 255;
      img(i, j, 1) = 255;
      img(i, j, 2) = 255;
    }
  }
  //img.display("edge", false);
  return img;
}

// 图像二值化处理
CImg<float> ImageSplit::convertToBinaryImg(CImg<float> warpingResult) {
  binaryImg = CImg<float>(warpingResult._width, warpingResult._height, 1, 1, 0);
  cimg_forXY(binaryImg, x, y) {
    int intensity = warpingResult(x, y, 0);
    //先去掉黑色边缘
    if (x <= BoundaryRemoveGap || y <= BoundaryRemoveGap
        || x >= warpingResult._width - BoundaryRemoveGap || y >= warpingResult._height - BoundaryRemoveGap) {
      binaryImg(x, y, 0) = 255; // 白色
    }
    else {
      if (intensity < BinaryGap)   // 小于阈值认为是边缘
        binaryImg(x, y, 0) = 0; 
      else
        binaryImg(x, y, 0) = 255;
    }
  }
  return binaryImg;
}

void ImageSplit::findDividingLine() {
  int lineColor[3]{255, 0, 0};
  histogramImg = CImg<float>(binaryImg._width, binaryImg._height, 1, 3, 255);
  dividingImg = binaryImg;
  vector<int> inflectionPoints; // 拐点
  cimg_forY(histogramImg, y) {
    int blackPixel = 0;
    // 统计黑点个数
    cimg_forX(binaryImg, x) {
      if (binaryImg(x, y, 0) == 0)
        blackPixel++;
      /*
      if (y < binaryImg._height-2 && binaryImg(x, y + 2, 0) == 0)
        blackPixel++;
        */
    }
    cimg_forX(histogramImg, x) {
      if (x < blackPixel) {
        histogramImg(x, y, 0) = 0;
        histogramImg(x, y, 1) = 0;
        histogramImg(x, y, 2) = 0;
      }
    }

    // 求Y方向直方图，谷的最少黑色像素个数为0
    // 判断是否为拐点
    if (y > 0) {
      // 下界
      if (blackPixel <= 0 && histogramImg(0, y - 1, 0) == 0) { 
        int yTemp = inflectionPoints[inflectionPoints.size() - 1];
        int count = 0;
        for (int i = yTemp; i < y; i++) {
          if (histogramImg(0, i, 0) == 0) {
            count++;
          }
        }
        
        // 统计两条线之间的黑点个数
        if (count >= 25 && (y - yTemp) > 1) {
          inflectionPoints.push_back(y);
         histogramImg.draw_line(0, y, histogramImg._width - 1, y, lineColor);
        }   
      }
      // 上界
      else if (blackPixel > 0 && histogramImg(0, y - 1, 0) != 0) {
        if (inflectionPoints.size() >= 1) {
          int yTemp = inflectionPoints[inflectionPoints.size() - 1];
          int count = 0;
          for (int i = yTemp; i < y; i++) {
            if (histogramImg(0, i, 0) == 0) {
              count++;
            }
          }
          // 统计两条线之间的黑点个数
          if (count == 0 && (y-1 - yTemp) > 2) {
            inflectionPoints.push_back(y - 1);
            histogramImg.draw_line(0, y - 1, histogramImg._width - 1, y - 1, lineColor);
          }
        }
        else {
          inflectionPoints.push_back(y - 1);
          histogramImg.draw_line(0, y - 1, histogramImg._width - 1, y - 1, lineColor);
        }
      }
    }
  }
  histogramImg.display("histogramImg Point");
  dividePoints.push_back(Point1(0, -1));
  // 两拐点中间做分割
  if (inflectionPoints.size() > 2) {
    for (int i = 1; i < inflectionPoints.size() - 1; i = i + 2) {
      int dividePoint = (inflectionPoints[i] + inflectionPoints[i + 1]) / 2;
      dividePoints.push_back(Point1(0, dividePoint));
      histogramImg.draw_line(0, dividePoint, histogramImg._width - 1, dividePoint, lineColor);
    }
  }
  //histogramImg.display("histogramImg Line");
  dividePoints.push_back(Point1(0, binaryImg._height - 1));
}

// 根据行分割线划分图片
void ImageSplit::divideIntoBarItemImg() {
  vector<Point1> tempDivideLinePointSet;
  for (int i = 1; i < dividePoints.size(); i++) {
    int barHeight = dividePoints[i].y - dividePoints[i - 1].y;
    int blackPixel = 0;
    // 初始化切割图
    CImg<float> barItemImg = CImg<float>(binaryImg._width, barHeight, 1, 1, 0);
    // 复制对应部分
    cimg_forXY(barItemImg, x, y) {
      barItemImg(x, y, 0) = binaryImg(x, dividePoints[i - 1].y + 1 + y, 0);
      if (barItemImg(x, y, 0) == 0)
        blackPixel++;
    }

    double blackPercent = (double)blackPixel / (double)(binaryImg._width * barHeight);
    // 只有当黑色像素个数超过图像大小一定比例0.005时，才可视作有数字
    if (blackPercent > 0.005) {
      //barItemImg.display(("barItemImg" + to_string(i)).c_str());
      string fileName = "./result/row/row";
      fileName = fileName + to_string(i) + ".bmp";
      barItemImg.save(fileName.c_str());

      // 横向切割
      vector<int> dividePosXset = getDivideLineXofSubImage(barItemImg);
      vector<CImg<float>> rowItemImgSet = getRowItemImgSet(barItemImg, dividePosXset);
      for (int j = 0; j < rowItemImgSet.size(); j++) {
        subImageSet.push_back(rowItemImgSet[j]);
        tempDivideLinePointSet.push_back(Point1(dividePosXset[j], dividePoints[i - 1].y));
      }
      if (i > 1) {
        histogramImg.draw_line(0, dividePoints[i - 1].y,
          histogramImg._width - 1, dividePoints[i - 1].y, lineColor);
        dividingImg.draw_line(0, dividePoints[i - 1].y,
          histogramImg._width - 1, dividePoints[i - 1].y, lineColor);
      }
      // 绘制竖线
      for (int j = 1; j < dividePosXset.size() - 1; j++) {
        dividingImg.draw_line(dividePosXset[j], dividePoints[i - 1].y,
          dividePosXset[j], dividePoints[i].y, lineColor);
      }
    }
  }
  //dividingImg.display("dividingImg");
  dividePoints.clear();
  for (int i = 0; i < tempDivideLinePointSet.size(); i++) {
    dividePoints.push_back(tempDivideLinePointSet[i]);
  }
}

// 膨胀
void ImageSplit::dilateImg(int barItemIndex) {
  // 膨胀Dilation -X-X-X-XYY方向
  CImg<float> answerXXY = CImg<float>(subImageSet[barItemIndex]._width,
                                        subImageSet[barItemIndex]._height, 1, 1, 0);
  cimg_forXY(subImageSet[barItemIndex], x, y) {
    answerXXY(x, y, 0) = getDilateXXY(subImageSet[barItemIndex], x, y);
  }

  // 膨胀Dilation -X-X-X-XYY方向
  CImg<float> answerXXY2 = CImg<float>(answerXXY._width, answerXXY._height, 1, 1, 0);
  cimg_forXY(answerXXY, x, y) {
    answerXXY2(x, y, 0) = getDilateXXY(answerXXY, x, y);
  }

  // 膨胀Dilation XY方向
  CImg<float> answerXY = CImg<float>(answerXXY2._width, answerXXY2._height, 1, 1, 0);
  cimg_forXY(answerXXY2, x, y) {
    answerXY(x, y, 0) = getDilateXY(answerXXY2, x, y);
  }

  cimg_forXY(subImageSet[barItemIndex], x, y) {
    subImageSet[barItemIndex](x, y, 0) = answerXY(x, y, 0);
  }
}

// 连通区域标记算法
void ImageSplit::connectedRegionsTagging(int barItemIndex) {
  tagImg = CImg<float>(subImageSet[barItemIndex]._width, subImageSet[barItemIndex]._height, 1, 1, 0);
  tagAccumulate = -1;

  cimg_forX(subImageSet[barItemIndex], x)
    cimg_forY(subImageSet[barItemIndex], y) {
    // 第一行和第一列
    if (x == 0 || y == 0) {
      int intensity = subImageSet[barItemIndex](x, y, 0);
      if (intensity == 0) {
        addNewTag(x, y, barItemIndex);
      }
    }
    // 其余的行和列
    else {
      int intensity = subImageSet[barItemIndex](x, y, 0);
      if (intensity == 0) {
        // 检查正上、左上、左中、左下这四个邻点
        int minTag = INT_MAX; //最小的tag
        Point1 minTagPointPos(-1, -1);
        // 先找最小的标记
        findMinTag(x, y, minTag, minTagPointPos, barItemIndex);

        // 当正上、左上、左中、左下这四个邻点有黑色点时，与minTag合并；
        if (minTagPointPos.x != -1 && minTagPointPos.y != -1) {
          mergeTagImageAndList(x, y - 1, minTag, minTagPointPos, barItemIndex);
          for (int i = -1; i <= 1; i++) {
            if (y + i < subImageSet[barItemIndex]._height)
              mergeTagImageAndList(x - 1, y + i, minTag, minTagPointPos, barItemIndex);
          }
          // 当前位置
          tagImg(x, y, 0) = minTag;
          Point1 cPoint(x + dividePoints[barItemIndex].x + 1, y + dividePoints[barItemIndex].y + 1);
          pointPosListSet[minTag].push_back(cPoint);

        }
        // 否则，作为新类
        else {
          addNewTag(x, y, barItemIndex);
        }
      }
    }
  }
}

// 存储分割后每一张数字的图以及对应的文件名称
void ImageSplit::saveSingleNumImg(int barItemIndex) {
  // 先统计每张数字图像黑色像素个数平均值
  int totalBlacks = 0, numberCount = 0;
  for (int i = 0; i < pointPosListSet.size(); i++) {
    if (pointPosListSet[i].size() != 0) {
      totalBlacks += pointPosListSet[i].size();
      numberCount++;
    }
  }
  int avgBlacks = totalBlacks / numberCount;

  int k = 0;
  for (int i = 0; i < pointPosListSet.size(); i++) {
    // 只有黑色像素个数大于平均值的一定比例，才可视为数字图像
    // 单张数字图像黑色像素个数超过所有数字图像，黑色像素个数均值的一定比例0.35才算作有数字
    if (pointPosListSet[i].size() != 0 && pointPosListSet[i].size() > avgBlacks * 0.35) {
      // 先找到数字的包围盒
      int xMin, xMax, yMin, yMax;
      getBoundingOfSingleNum(i, xMin, xMax, yMin, yMax);

      int width = xMax - xMin;
      int height = yMax - yMin;

      // 将单个数字填充到新图像：扩充到正方形
      // 单张数字图像边缘填充宽度为5
      int imgSize = (width > height ? width : height) + 5 * 2;
      CImg<float> singleNum = CImg<float>(imgSize, imgSize, 1, 1, 0);

      list<Point1>::iterator it = pointPosListSet[i].begin();
      for (; it != pointPosListSet[i].end(); it++) {
        int x = (*it).x;
        int y = (*it).y;
        int singleNumImgPosX, singleNumImgPosY;
        if (height > width) {
          singleNumImgPosX = (x - xMin) + (imgSize - width) / 2;
          singleNumImgPosY = (y - yMin) + 5;
        }
        else {
          singleNumImgPosX = (x - xMin) + 5;
          singleNumImgPosY = (y - yMin) + (imgSize - height) / 2;
        }
        singleNum(singleNumImgPosX, singleNumImgPosY, 0) = 255;
      }

      // 对单个数字图像做Y方向腐蚀操作
      singleNum = eroseImg(singleNum);

      string postfix = ".bmp";
      char shortImgName[200];
      //sprintf(shortImgName, "%d_%d%s\n", barItemIndex, classTagSet[i], postfix.c_str());
      sprintf(shortImgName, "%d_%d%s\n", barItemIndex, k, postfix.c_str());
      imglisttxt += string(shortImgName);

      char addr[200];
      //sprintf(addr, "%s%d_%d%s", basePath.c_str(), barItemIndex, classTagSet[i], postfix.c_str());
      sprintf(addr, "%s%d_%d%s", basePath.c_str(), barItemIndex, k, postfix.c_str());
      singleNum.save(addr);
      k++;
      pointPosListSetForDisplay.push_back(pointPosListSet[i]);
    }
  }
  // 存储图像文件名
  imglisttxt += "*\n";
  ofstream predictImageListOutput(basePath + "imagelist.txt");
  predictImageListOutput << imglisttxt.c_str();
  predictImageListOutput.close();

  // 把tag集、每一类链表数据集清空
  classTagSet.clear();
  for (int i = 0; i < pointPosListSet.size(); i++) {
    pointPosListSet[i].clear();
  }
  pointPosListSet.clear();
}

// 获取每一行的每个数字的竖直分割线
vector<int> ImageSplit::getDivideLineXofSubImage(const CImg<float>& subImg) {
  // 先绘制X方向灰度直方图
  CImg<float> XHistogramImage = CImg<float>(subImg._width, subImg._height, 1, 3, 255);
  cimg_forX(subImg, x) {
    int blackPixel = 0;
    cimg_forY(subImg, y) {
      if (subImg(x, y, 0) == 0)
        blackPixel++;
    }
    // 对于每一列x，只有黑色像素多于一定值，才绘制在直方图上（确认是数字）
    // 求X方向直方图，谷的最少黑色像素个数        
    if (blackPixel >= 4) {
      cimg_forY(subImg, y) {
        if (y < blackPixel) {
          XHistogramImage(x, y, 0) = 0;
          XHistogramImage(x, y, 1) = 0;
          XHistogramImage(x, y, 2) = 0;
        }
      }
    }
  }

  vector<int> InflectionPosXs = getInflectionPosXs(XHistogramImage);    //获取拐点
  for (int i = 0; i < InflectionPosXs.size(); i++) {
    XHistogramImage.draw_line(InflectionPosXs[i], 0, InflectionPosXs[i], XHistogramImage._height - 1, lineColor);
  }

  // 两拐点中间做分割
  vector<int> dividePosXs;
  dividePosXs.push_back(-1);
  if (InflectionPosXs.size() > 2) {
    for (int i = 1; i < InflectionPosXs.size() - 1; i = i + 2) {
      int divideLinePointX = (InflectionPosXs[i] + InflectionPosXs[i + 1]) / 2;
      dividePosXs.push_back(divideLinePointX);
    }
  }
  dividePosXs.push_back(XHistogramImage._width - 1);
  return dividePosXs;
}

// 根据X方向直方图判断真实的拐点
vector<int> ImageSplit::getInflectionPosXs(const CImg<float>& XHistogramImage) {
  vector<int> resultInflectionPosXs, tempInflectionPosXs;
  int totalDist = 0, dist = 0;
  // 查找拐点
  //XHistogramImage.display("XHistogramImage");
  // 置反，与MNIST数据一致
  cimg_forX(XHistogramImage, x) {
    if (x >= 1) {
      // 白转黑
      if (XHistogramImage(x, 0, 0) == 0 && XHistogramImage(x - 1, 0, 0) == 255) 
        tempInflectionPosXs.push_back(x - 1);
      // 黑转白
      else if (XHistogramImage(x, 0, 0) == 255 && XHistogramImage(x - 1, 0, 0) == 0)
        tempInflectionPosXs.push_back(x);
    }
  }
  for (int i = 2; i < tempInflectionPosXs.size() - 1; i = i + 2) {
    int tempdist = tempInflectionPosXs[i] - tempInflectionPosXs[i - 1];
    if (tempdist <= 0)
      tempdist--;
    totalDist += tempdist;
  }

  // 计算间距平均距离
  dist += (tempInflectionPosXs.size() - 2) / 2;
  int avgDist = 0;
  if (dist != 0) {
    avgDist = totalDist / dist;
  }

  resultInflectionPosXs.push_back(tempInflectionPosXs[0]); //头
  // 当某个间距大于平均距离的一定倍数时，视为分割点所在间距
  for (int i = 2; i < tempInflectionPosXs.size() - 1; i = i + 2) {
    int dist = tempInflectionPosXs[i] - tempInflectionPosXs[i - 1];
    if (dist > avgDist * 4) {
      resultInflectionPosXs.push_back(tempInflectionPosXs[i - 1]);
      resultInflectionPosXs.push_back(tempInflectionPosXs[i]);
    }
  }
  resultInflectionPosXs.push_back(tempInflectionPosXs[tempInflectionPosXs.size() - 1]); //尾
  return resultInflectionPosXs;
}

// 分割行子图，得到单个数字图
vector<CImg<float>> ImageSplit::getRowItemImgSet(const CImg<float>& lineImg, vector<int> _dividePosXset) {
  vector<CImg<float>> result;
  
  for (int i = 1; i < _dividePosXset.size(); i++) {
    int rowItemWidth = _dividePosXset[i] - _dividePosXset[i - 1];
    // 初始化单个数字图
    CImg<float> rowItemImg = CImg<float>(rowItemWidth, lineImg._height, 1, 1, 0);
    // 复制对应部分
    cimg_forXY(rowItemImg, x, y) {
      rowItemImg(x, y, 0) = lineImg(x + _dividePosXset[i - 1] + 1, y, 0);
    }
    result.push_back(rowItemImg);
    //rowItemImg.display(("rowItemImg" + to_string(i)).c_str());
    
  }
  return result;
}

// XY方向的正扩张
int ImageSplit::getDilateXY(const CImg<float>& Img, int x, int y) {
  int intensity = Img(x, y, 0);
  if (intensity == 255) { // 若中间点为白色
    // X和Y都是一个单位
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        // 保证不越界
        if (0 <= x + i && x + i < Img._width && 0 <= y + j && y + j < Img._height) {
          // 水平或垂直方向
          if (i != -1 && j != -1 || i != 1 && j != 1 || i != 1 && j != -1 || i != -1 && j != 1)
            if (Img(x + i, y + j, 0) == 0) {
              intensity = 0; // 扩张
              break;
            }
        }
      }
      if (intensity != 255) {
        break;
      }
    }
  }
  return intensity;
}

// X方向2个单位的负扩张，Y方向1个单位的正扩张
int ImageSplit::getDilateXXY(const CImg<float>& Img, int x, int y) {
  int intensity = Img(x, y, 0);
  if (intensity == 255) {    //若中间点为白色
    int blackAccu = 0;
    // 一个单位
    for (int i = -1; i <= 1; i++) {
      if (0 <= y + i && y + i < Img._height) {    //竖直方向
        if (Img(x, y + i, 0) == 0)
          blackAccu++;
      }
    }
    // 两个单位
    for (int i = -2; i <= 2; i++) {
      if (0 <= x + i && x + i < Img._width) {     //水平方向
        if (Img(x + i, y, 0) == 0)
          blackAccu--;
      }
    }
    if (blackAccu > 0) {
      intensity = 0;  // 扩张
    }
  }
  return intensity;
}

// 添加新的类tag
void ImageSplit::addNewTag(int x, int y, int barItemIndex) {
  tagAccumulate++;
  tagImg(x, y, 0) = tagAccumulate;
  classTagSet.push_back(tagAccumulate);
  list<Point1> pList;
  Point1 cPoint(x + dividePoints[barItemIndex].x + 1, y + dividePoints[barItemIndex].y + 1);
  pList.push_back(cPoint);
  pointPosListSet.push_back(pList);
}

// 在正上、左上、正左、左下这四个邻点中找到最小的tag
void ImageSplit::findMinTag(int x, int y, int &minTag, Point1 &minTagPointPos, int barItemIndex) {
  // 正上
  if (subImageSet[barItemIndex](x, y - 1, 0) == 0) {
    if (tagImg(x, y - 1, 0) < minTag) {
      minTag = tagImg(x, y - 1, 0);
      minTagPointPos.x = x;
      minTagPointPos.y = y - 1;
    }
  }
  // 左上、左中、左下
  for (int i = -1; i <= 1; i++) {
    if (y + i < subImageSet[barItemIndex]._height) {
      if (subImageSet[barItemIndex](x - 1, y + i, 0) == 0 && tagImg(x - 1, y + i, 0) < minTag) {
        minTag = tagImg(x - 1, y + i, 0);
        minTagPointPos.x = x - 1;
        minTagPointPos.y = y + i;
      }
    }
  }
}

// 合并某个点(x,y)所属类别
void ImageSplit::mergeTagImageAndList(int x, int y, const int minTag, const Point1 minTagPointPos, int barItemIndex) {
  // 赋予最小标记，合并列表
  if (subImageSet[barItemIndex](x, y, 0) == 0) {
    int tagBefore = tagImg(x, y, 0);
    if (tagBefore != minTag) {
      //把所有同一类的tag替换为最小tag、把list接到最小tag的list
      list<Point1>::iterator it = pointPosListSet[tagBefore].begin();
      for (; it != pointPosListSet[tagBefore].end(); it++) {
        tagImg((*it).x - dividePoints[barItemIndex].x - 1, (*it).y - dividePoints[barItemIndex].y - 1, 0) = minTag;
      }
      pointPosListSet[minTag].splice(pointPosListSet[minTag].end(), pointPosListSet[tagBefore]);
    }
  }
}

// 获取单个数字的包围盒
void ImageSplit::getBoundingOfSingleNum(int listIndex, int& xMin, int& xMax, int& yMin, int& yMax) {
  xMin = yMin = INT_MAX;
  xMax = yMax = -1;
  if (!pointPosListSet.empty()) {
    list<Point1>::iterator it = pointPosListSet[listIndex].begin();
    for (; it != pointPosListSet[listIndex].end(); it++) {
      int x = (*it).x, y = (*it).y;
      xMin = x < xMin ? x : xMin;
      yMin = y < yMin ? y : yMin;
      xMax = x > xMax ? x : xMax;
      yMax = y > yMax ? y : yMax;
    }
  }
  else {
    list<Point1>::iterator it = pointPosListSetForDisplay[listIndex].begin();
    for (; it != pointPosListSetForDisplay[listIndex].end(); it++) {
      int x = (*it).x, y = (*it).y;
      xMin = x < xMin ? x : xMin;
      yMin = y < yMin ? y : yMin;
      xMax = x > xMax ? x : xMax;
      yMax = y > yMax ? y : yMax;
    }
  }
}

// 对单个数字图像做Y方向腐蚀操作
CImg<float> ImageSplit::eroseImg(CImg<float>& Img) {
  CImg<float> result = CImg<float>(Img._width, Img._height, 1, 1, 0);
  cimg_forXY(Img, x, y) {
    result(x, y, 0) = Img(x, y, 0);
    if (Img(x, y, 0) == 255) {
      // 上方
      if (y - 1 >= 0) {
        if (Img(x, y - 1, 0) == 0)
          result(x, y, 0) = 0;
      }
      // 下方
      if (y + 1 < Img._height) {
        if (Img(x, y + 1, 0) == 0)
          result(x, y, 0) = 0;
      }
    }
  }
  return result;
}