
#include "ImageWrap.h"
#include "Blend.h"
#include "Stitch.h"


// 获得灰度图，加快处理速度
CImg<float> get_gray_image(CImg<float>& image) {
  CImg<float> res(image._width, image._height);
  cimg_forXY(image, x, y) {
    res(x, y) = 0.299 * image(x, y, 0, 0) +
                0.587 * image(x, y, 0, 1) +
                0.114 * image(x, y, 0, 2);
  }
  return res;
}

// 将A转到B，且x,y交换（即交换拼接方向）
void ReplacePairs(vector<POINT_PAIR>& A, vector<POINT_PAIR>& B) {
  B.clear();
  for (int i = 0; i < A.size(); i++) {
    // src 与 dst交换
    POINT_PAIR temp(A[i].b, A[i].a);
    B.push_back(temp);
  }
}

// 获取平均位移
vector<int> getAvgOffset(const vector<POINT_PAIR>& pairs, vector<int>& indices) {
  int offset_x = 0;
  int offset_y = 0;
  int min_x = 1000;
  int min_y = 1000;

  int size = indices.size();
  int cnt = 0;
  for (int i = 0; i < size; i++) {

    int diff_x = pairs[i].a.x - pairs[i].b.x;
    int diff_y = pairs[i].a.y - pairs[i].b.y;
    if (diff_x == 0 || diff_y == 0 || abs(diff_y) > 220) {
      continue;
    }
    // 求出位移的和
    offset_x += diff_x;
    offset_y += diff_y;
    // 统计点的数量
    cnt++;

    // 最小的x
    if (pairs[i].a.x < min_x) {
      min_x = pairs[i].a.x;
    }
    
    // 最小的y
    if (pairs[i].a.y < min_y) {
      min_y = pairs[i].a.y;
    }
  }

  // 求平均值
  offset_x /= cnt;
  offset_y /= cnt;

  int ans_x = 0;
  int ans_y = 0;
  cnt = 0;
  for (int i = 0; i < size; i++) {
    int diff_x = pairs[i].a.x - pairs[i].b.x;
    int diff_y = pairs[i].a.y - pairs[i].b.y;
    cout << "diff_x " << diff_x << " diff_y " << diff_y << endl;
    

    if (abs(diff_x - offset_x) < abs(offset_x) / 2 && (abs(diff_y - offset_y) < abs(offset_y) / 2 || abs(offset_y) / 2 < 4)) {
      ans_x += diff_x;
      ans_y += diff_y;
      cnt++;
    }
  }

  if (cnt) {
    ans_x /= cnt;
    ans_y /= cnt;
  }
  vector<int> res;
  res.push_back(ans_x);
  res.push_back(ans_y);
  res.push_back(min_x);
  res.push_back(min_y);
  return res;
}

// 使用SIFT提取特征点
map<vector<float>, VlSiftKeypoint> extractFeatures(CImg<float>& img) {
  CImg<float> src(img);
  float resize_factor;
  int width = src._width;
  int height = src._height;

  // 优化计算速度
  if (width < height) {
    resize_factor = RESIZE_SIZE / width;
  }
  else {
    resize_factor = RESIZE_SIZE / height;
  }

  if (resize_factor >= 1) {
    resize_factor = 1;
  }
  else {
    src.resize(width * resize_factor, height * resize_factor, src.spectrum(), 3);
  }

  // vl_sift_pix 就是float型数据
  vl_sift_pix *imageData = new vl_sift_pix[src._height * src._width];

  // 设置SIFT算法过滤参数
  // Setting SIFT filter params
  int noctaves = 4, nlevels = 2, o_min = 0;

  // 图像二维转一维
  for (int i = 0; i < src.width(); i++) {
    for (int j = 0; j < src.height(); j++) {
      imageData[j * src.width() + i] = src(i, j, 0);
    }
  }

  // 这个过滤器实现了SIFT检测器和描述符
  VlSiftFilt *sf = NULL;

  // noctaves: numbers of octaves 组数
  // nlevels: numbers of levels per octave 每组的层数
  // o_min: first octave index 第一组的索引号

  sf = vl_sift_new(src.width(), src.height(), noctaves, nlevels, o_min);

  map<vector<float>, VlSiftKeypoint> features; // 记录特征点的描述符，一个特征点有可能有多个描述符，最多有4个

                                               // Compute the first octave of the DOG scale space
                                               // 这个函数开始处理一幅新图像，通过计算它在低层的高斯尺度空间
                                               // 它还清空内部记录关键点的缓冲区
  if (vl_sift_process_first_octave(sf, imageData) != VL_ERR_EOF) {
    while (1) {
      // Run the SIFT detector to get the keypoints.
      vl_sift_detect(sf); // 计算每组中的关键点

      VlSiftKeypoint *pKeyPoint = sf->keys;

      // 遍历每个特征点
      for (int i = 0; i < sf->nkeys; i++) {
        VlSiftKeypoint tempKp = *pKeyPoint;

        // 计算并遍历每个点的方向
        double angles[4];

        // 计算每个极值点的方向，包括主方向和辅方向，最多4个方向
        int angleCount = vl_sift_calc_keypoint_orientations(sf, angles, &tempKp); // 方向数量

        for (int j = 0; j < angleCount; j++) {
          // 计算每个方向的描述符
          vl_sift_pix descriptors[128];

          // 获取特征点的描述子
          vl_sift_calc_keypoint_descriptor(sf, descriptors, &tempKp, angles[j]);

          // 复制到vector
          vector<float> des;
          int k = 0;
          while (k < 128) {
            des.push_back(descriptors[k]);
            k++;
          }

          // 处理特征点信息
          tempKp.x /= resize_factor;
          tempKp.y /= resize_factor;
          tempKp.ix = tempKp.x;
          tempKp.iy = tempKp.y;

          features.insert(make_pair(des, tempKp)); // 插入到特征点map
        }

        pKeyPoint++;
      }
      // 这个函数计算高斯尺度空间中的下一组尺度空间图像
      // 这个函数会清除在前一层空间中检测到的特征点
      if (vl_sift_process_next_octave(sf) == VL_ERR_EOF) {
        break;
      }
    }
  }

  // 释放资源
  vl_sift_delete(sf);
  delete[] imageData;
  imageData = NULL;

  return features;
}


// 缝合一系列图片
CImg<float> stitching(vector<CImg<float> >& src_imgs) {
  // 存储特征值等数据
  vector<map<vector<float>, VlSiftKeypoint> > features(src_imgs.size());

  // 对原素材进行预处理
  for (int i = 0; i < src_imgs.size(); i++) {
    // 柱面投影
    cout << "Cylinder Projection" << endl;
    src_imgs[i] = CylinderProjection(src_imgs[i]);

    // 转化为灰度图
    CImg<float> gray = get_gray_image(src_imgs[i]);

    // 对每个灰度图进行特征提取
    cout << "SIFT: Extract Features" << endl;
    features[i] = extractFeatures(gray);
  }

  // 邻居表
  bool adjacent[20][20] = { false };
  vector<vector<int> > matching_index(src_imgs.size());


  // 找到每张图片的邻居，确认缝合对象
  cout << "Find Adjacent images.\n";
  for (int i = 0; i < src_imgs.size(); i++) {
    for (int j = i + 1; j < src_imgs.size(); j++) {
      // 对比两幅图的特征点，求出match的特征点集合
      vector<POINT_PAIR> pairs = getPointPairsFromFeatures(features[i], features[j]);

      // 如果吻合点数量超过30，则认为两幅图是相邻的
      if (pairs.size() >= 20) {
        // 记录相邻关系
        adjacent[i][j] = true;
        matching_index[i].push_back(j);
        cout << "Adjacent: " << i << " and " << j << endl;
      }
    }
  }
  cout << endl;

  cout << "Stitching" << endl;
  int beginIndex = 0;

  // 待拼接队列
  queue<int> unstitched_idx;
  unstitched_idx.push(beginIndex);

  // 当前已拼接图片
  CImg<float> cur_stitched_img = src_imgs[beginIndex];

  while (!unstitched_idx.empty()) {
    int sourceIndex = unstitched_idx.front();
    unstitched_idx.pop();

    for (int i = 0; i < matching_index[sourceIndex].size(); i--) {
      // 与当前图片拼接的图片下标
      int nextIndex = matching_index[sourceIndex][i];

      
      if (adjacent[sourceIndex][nextIndex]) {
        adjacent[sourceIndex][nextIndex] = adjacent[nextIndex][sourceIndex] = false;
        unstitched_idx.push(nextIndex);


        cout << "get Features.\n";
        // kd树找最近邻
        vector<POINT_PAIR> src_to_dst_pairs = getPointPairsFromFeatures(features[sourceIndex], features[nextIndex]);
        vector<POINT_PAIR> dst_to_src_pairs = getPointPairsFromFeatures(features[nextIndex], features[sourceIndex]);
        // 找最佳匹配方向
        if (src_to_dst_pairs.size() > dst_to_src_pairs.size())
          ReplacePairs(src_to_dst_pairs, dst_to_src_pairs);
        else
          ReplacePairs(dst_to_src_pairs, src_to_dst_pairs);

        // RANSAC算法
        cout << "RANSAC" << endl;
        vector<int> indices = RANSAC(dst_to_src_pairs);

        // 根据最佳模型，计算平均位移
        vector<int> offset = getAvgOffset(src_to_dst_pairs, indices);

        cout << "offset_x " << offset[0] << " offset_y" << offset[1] << endl;

        // 使用平均位移信息拼接图片
        cur_stitched_img = blending(cur_stitched_img, src_imgs[nextIndex], offset[0], offset[1], offset[2], offset[3]);
        cur_stitched_img.display("mid-process", false);
      }
    }
  }
  return cur_stitched_img;
}