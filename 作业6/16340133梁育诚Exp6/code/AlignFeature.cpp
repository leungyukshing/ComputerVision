#include "AlignFeature.h"
#include "ImageWrap.h"

// 计算迭代次数
int numberOfIteration(float p, float w, int n) {
  return ceil(log(1 - p) / log(1 - pow(w, n)));
}

// 产生min-max之间的随机数
int random(int min, int max) {
  return rand() % (max - min + 1) + min;
}

// 从两幅图像的特征点集合中找出匹配的点对
vector<POINT_PAIR> getPointPairsFromFeatures(map<vector<float>, VlSiftKeypoint>& feature_a,
    map<vector<float>, VlSiftKeypoint>& feature_b) {
  // 使用KD树来寻找匹配点集
  VlKDForest* forest = vl_kdforest_new(VL_TYPE_FLOAT, 128, 1, VlDistanceL1);

  // 提取出128维的特征向量
  float* data = new float[128 * feature_a.size()];
  int k = 0;
  for (auto it = feature_a.begin(); it != feature_a.end(); it++) {
    const vector<float>& descriptors = it->first;

    for (int i = 0; i < 128; i++) {
      data[i + 128 * k] = descriptors[i];
    }
    k++;
  }

  // 构建kd树
  vl_kdforest_build(forest, feature_a.size(), data);

  vector<POINT_PAIR> result;
  // 构建一个searcher
  VlKDForestSearcher* searcher = vl_kdforest_new_searcher(forest);
  VlKDForestNeighbor neighbours[2];

  for (auto it = feature_b.begin(); it != feature_b.end(); it++) {
    float *temp_data = new float[128];

    for (int i = 0; i < 128; i++) {
      temp_data[i] = (it->first)[i];
    }

    // 找最近的两个近邻的距离
    int nvisited = vl_kdforestsearcher_query(searcher, neighbours, 2, temp_data);

    // 两个邻居的距离比值（达朗贝尔判别法）
    float ratio = neighbours[0].distance / neighbours[1].distance;
    if (ratio < 0.6) {
      vector<float> des(128);
      for (int j = 0; j < 128; j++) {
        des[j] = data[j + neighbours[0].index * 128];
      }

      VlSiftKeypoint left = feature_a.find(des)->second;
      VlSiftKeypoint right = it->second;
      result.push_back(POINT_PAIR(left, right));
    }

    delete[] temp_data;
    temp_data = NULL;
  }
  vl_kdforestsearcher_delete(searcher);
  vl_kdforest_delete(forest);

  delete[] data;
  data = NULL;

  return result;
}


vector<int> getIndicesOfInlier(const vector<POINT_PAIR>& pairs, MatrixXf& H, set<int>& selected_indices) {
  vector<int> inliner_indices;

  for (int i = 0; i < pairs.size(); i++) {
    // 跳过已选择的点
    if (selected_indices.find(i) != selected_indices.end())
      continue;

    float real_x = pairs[i].b.x;
    float real_y = pairs[i].b.y;

    if (real_x == 0 || real_y == 0)
      continue;

    // 根据透视变换矩阵得到点坐标，验证模型
    float x = getXAfterWarping(pairs[i].a.x, pairs[i].a.y, H);
    float y = getYAfterWarping(pairs[i].a.x, pairs[i].a.y, H);

    // 小于一定距离的视作局内点
    float dist = sqrt((x - real_x) * (x - real_x) + (y - real_y) * (y - real_y));
    if (dist < RANSAC_THRESHOLD) {
      inliner_indices.push_back(i);
    }
  }
  return inliner_indices;
}

// 根据四对点，计算透视变换矩阵
MatrixXf getHomographyFromPointPairs(vector<POINT_PAIR>& pairs) {
  if (pairs.size() != 4) {
    cout << "Fail to get Homography matrix.n";
    exit(1);
  }

  MatrixXf UV(8, 1), A(8, 8);
  UV << pairs[0].b.x, pairs[0].b.y,
  pairs[1].b.x, pairs[1].b.y,
  pairs[2].b.x, pairs[2].b.y,
  pairs[3].b.x, pairs[3].b.y;

  A << pairs[0].a.x, pairs[0].a.y, 1, 0, 0, 0, -pairs[0].a.x * pairs[0].b.x, -pairs[0].a.y * pairs[0].b.x,
  0, 0, 0, pairs[0].a.x, pairs[0].a.y, 1, -pairs[0].a.x * pairs[0].b.y, -pairs[0].a.y * pairs[0].b.y,
  pairs[1].a.x, pairs[1].a.y, 1, 0, 0, 0, -pairs[1].a.x * pairs[1].b.x, -pairs[1].a.y * pairs[1].b.x,
  0, 0, 0, pairs[1].a.x, pairs[1].a.y, 1, -pairs[1].a.x * pairs[1].b.y, -pairs[1].a.y * pairs[1].b.y,
  pairs[2].a.x, pairs[2].a.y, 1, 0, 0, 0, -pairs[2].a.x * pairs[2].b.x, -pairs[2].a.y * pairs[2].b.x,
  0, 0, 0, pairs[2].a.x, pairs[2].a.y, 1, -pairs[2].a.x * pairs[2].b.y, -pairs[2].a.y * pairs[2].b.y,
  pairs[3].a.x, pairs[3].a.y, 1, 0, 0, 0, -pairs[3].a.x * pairs[3].b.x, -pairs[3].a.y * pairs[3].b.x,
  0, 0, 0, pairs[3].a.x, pairs[3].a.y, 1, -pairs[3].a.x * pairs[3].b.y, -pairs[3].a.y * pairs[3].b.y;

  A = A.inverse();
  // 透视变换矩阵
  MatrixXf M(8, 1);
  M = A * UV;
  return M;
}
/*
// 最小二乘对局内点拟合
MatrixXf leastSquareSolution(const vector<POINT_PAIR>& pairs, vector<int>& idxs) {
  int calc_size = idxs.size();

  CImg<double> A(4, calc_size, 1, 1, 0);
  CImg<double> b(1, calc_size, 1, 1, 0);

  for (int i = 0; i < calc_size; i++) {
    int cur_index = idxs[i];

    A(0, i) = pairs[cur_index].a.x;
    A(1, i) = pairs[cur_index].a.y;
    A(2, i) = pairs[cur_index].a.x * pairs[cur_index].a.y;
    A(3, i) = 1;

    b(0, i) = pairs[cur_index].b.x;
  }

  CImg<double> x1 = b.get_solve(A);

  for (int i = 0; i < calc_size; i++) {
    int cur_index = idxs[i];

    b(0, i) = pairs[cur_index].b.y;
  }

  CImg<double> x2 = b.get_solve(A);

  MatrixXf H(8, 1);
  H << x1(0, 0), x1(0, 1), x1(0, 2), x1(0, 3), x2(0, 0), x2(0, 1), x2(0, 2), x2(0, 3);

  return H;
}
*/

// RANSAC算法
vector<int> RANSAC(const vector<POINT_PAIR>& pairs) {
  if (pairs.size() < NUM_OF_PAIR) {
    cout << "Not enough pairs.\n";
    exit(1);
  }

  // 初始化随机数种子
  srand(time(0));
  int iterations = numberOfIteration(CONFIDENCE, INLINER_RATIO, NUM_OF_PAIR);

  // 最佳模型的下标
  vector<int> max_inliner_indices;

  // 抽样迭代
  while (iterations--) {
    vector<POINT_PAIR> random_point_pairs; // 随机点集
    // 每次选择的集合
    set<int> selected_indices;

    // 1. 随机选取四个点
    for (int i = 0; i < NUM_OF_PAIR; i++) {
      int idx = random(0, pairs.size() - 1);
      // All pairs are unique
      while (selected_indices.find(idx) != selected_indices.end()) {
        idx = random(0, pairs.size() - 1);
      }
      selected_indices.insert(idx); // 存入已选择下集合

      random_point_pairs.push_back(pairs[idx]); // 将该点放入随即点集
    }

    // 2. 求解变换模型
    MatrixXf H = getHomographyFromPointPairs(random_point_pairs);

    // 3. 接受局内点，拒绝局外点，SSD(pi’, H pi)< ε
    vector<int> cur_inliner_indices = getIndicesOfInlier(pairs, H, selected_indices);

    // 4. 保留投票数最多的模型（含局内点最多的H）
    if (cur_inliner_indices.size() > max_inliner_indices.size())
      max_inliner_indices = cur_inliner_indices;
  }
  // 递归结束
  // 5. 用所有假设的局内点去重新估计模型（使用最小二乘法）
  //MatrixXf H = leastSquareSolution(pairs, max_inliner_indices);

  return max_inliner_indices;
}

