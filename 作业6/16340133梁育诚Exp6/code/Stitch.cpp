
#include "ImageWrap.h"
#include "Blend.h"
#include "Stitch.h"


// ��ûҶ�ͼ���ӿ촦���ٶ�
CImg<float> get_gray_image(CImg<float>& image) {
  CImg<float> res(image._width, image._height);
  cimg_forXY(image, x, y) {
    res(x, y) = 0.299 * image(x, y, 0, 0) +
                0.587 * image(x, y, 0, 1) +
                0.114 * image(x, y, 0, 2);
  }
  return res;
}

// ��Aת��B����x,y������������ƴ�ӷ���
void ReplacePairs(vector<POINT_PAIR>& A, vector<POINT_PAIR>& B) {
  B.clear();
  for (int i = 0; i < A.size(); i++) {
    // src �� dst����
    POINT_PAIR temp(A[i].b, A[i].a);
    B.push_back(temp);
  }
}

// ��ȡƽ��λ��
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
    // ���λ�Ƶĺ�
    offset_x += diff_x;
    offset_y += diff_y;
    // ͳ�Ƶ������
    cnt++;

    // ��С��x
    if (pairs[i].a.x < min_x) {
      min_x = pairs[i].a.x;
    }
    
    // ��С��y
    if (pairs[i].a.y < min_y) {
      min_y = pairs[i].a.y;
    }
  }

  // ��ƽ��ֵ
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

// ʹ��SIFT��ȡ������
map<vector<float>, VlSiftKeypoint> extractFeatures(CImg<float>& img) {
  CImg<float> src(img);
  float resize_factor;
  int width = src._width;
  int height = src._height;

  // �Ż������ٶ�
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

  // vl_sift_pix ����float������
  vl_sift_pix *imageData = new vl_sift_pix[src._height * src._width];

  // ����SIFT�㷨���˲���
  // Setting SIFT filter params
  int noctaves = 4, nlevels = 2, o_min = 0;

  // ͼ���άתһά
  for (int i = 0; i < src.width(); i++) {
    for (int j = 0; j < src.height(); j++) {
      imageData[j * src.width() + i] = src(i, j, 0);
    }
  }

  // ���������ʵ����SIFT�������������
  VlSiftFilt *sf = NULL;

  // noctaves: numbers of octaves ����
  // nlevels: numbers of levels per octave ÿ��Ĳ���
  // o_min: first octave index ��һ���������

  sf = vl_sift_new(src.width(), src.height(), noctaves, nlevels, o_min);

  map<vector<float>, VlSiftKeypoint> features; // ��¼���������������һ���������п����ж���������������4��

                                               // Compute the first octave of the DOG scale space
                                               // ���������ʼ����һ����ͼ��ͨ���������ڵͲ�ĸ�˹�߶ȿռ�
                                               // ��������ڲ���¼�ؼ���Ļ�����
  if (vl_sift_process_first_octave(sf, imageData) != VL_ERR_EOF) {
    while (1) {
      // Run the SIFT detector to get the keypoints.
      vl_sift_detect(sf); // ����ÿ���еĹؼ���

      VlSiftKeypoint *pKeyPoint = sf->keys;

      // ����ÿ��������
      for (int i = 0; i < sf->nkeys; i++) {
        VlSiftKeypoint tempKp = *pKeyPoint;

        // ���㲢����ÿ����ķ���
        double angles[4];

        // ����ÿ����ֵ��ķ��򣬰���������͸��������4������
        int angleCount = vl_sift_calc_keypoint_orientations(sf, angles, &tempKp); // ��������

        for (int j = 0; j < angleCount; j++) {
          // ����ÿ�������������
          vl_sift_pix descriptors[128];

          // ��ȡ�������������
          vl_sift_calc_keypoint_descriptor(sf, descriptors, &tempKp, angles[j]);

          // ���Ƶ�vector
          vector<float> des;
          int k = 0;
          while (k < 128) {
            des.push_back(descriptors[k]);
            k++;
          }

          // ������������Ϣ
          tempKp.x /= resize_factor;
          tempKp.y /= resize_factor;
          tempKp.ix = tempKp.x;
          tempKp.iy = tempKp.y;

          features.insert(make_pair(des, tempKp)); // ���뵽������map
        }

        pKeyPoint++;
      }
      // ������������˹�߶ȿռ��е���һ��߶ȿռ�ͼ��
      // ��������������ǰһ��ռ��м�⵽��������
      if (vl_sift_process_next_octave(sf) == VL_ERR_EOF) {
        break;
      }
    }
  }

  // �ͷ���Դ
  vl_sift_delete(sf);
  delete[] imageData;
  imageData = NULL;

  return features;
}


// ���һϵ��ͼƬ
CImg<float> stitching(vector<CImg<float> >& src_imgs) {
  // �洢����ֵ������
  vector<map<vector<float>, VlSiftKeypoint> > features(src_imgs.size());

  // ��ԭ�زĽ���Ԥ����
  for (int i = 0; i < src_imgs.size(); i++) {
    // ����ͶӰ
    cout << "Cylinder Projection" << endl;
    src_imgs[i] = CylinderProjection(src_imgs[i]);

    // ת��Ϊ�Ҷ�ͼ
    CImg<float> gray = get_gray_image(src_imgs[i]);

    // ��ÿ���Ҷ�ͼ����������ȡ
    cout << "SIFT: Extract Features" << endl;
    features[i] = extractFeatures(gray);
  }

  // �ھӱ�
  bool adjacent[20][20] = { false };
  vector<vector<int> > matching_index(src_imgs.size());


  // �ҵ�ÿ��ͼƬ���ھӣ�ȷ�Ϸ�϶���
  cout << "Find Adjacent images.\n";
  for (int i = 0; i < src_imgs.size(); i++) {
    for (int j = i + 1; j < src_imgs.size(); j++) {
      // �Ա�����ͼ�������㣬���match�������㼯��
      vector<POINT_PAIR> pairs = getPointPairsFromFeatures(features[i], features[j]);

      // ����Ǻϵ���������30������Ϊ����ͼ�����ڵ�
      if (pairs.size() >= 20) {
        // ��¼���ڹ�ϵ
        adjacent[i][j] = true;
        matching_index[i].push_back(j);
        cout << "Adjacent: " << i << " and " << j << endl;
      }
    }
  }
  cout << endl;

  cout << "Stitching" << endl;
  int beginIndex = 0;

  // ��ƴ�Ӷ���
  queue<int> unstitched_idx;
  unstitched_idx.push(beginIndex);

  // ��ǰ��ƴ��ͼƬ
  CImg<float> cur_stitched_img = src_imgs[beginIndex];

  while (!unstitched_idx.empty()) {
    int sourceIndex = unstitched_idx.front();
    unstitched_idx.pop();

    for (int i = 0; i < matching_index[sourceIndex].size(); i--) {
      // �뵱ǰͼƬƴ�ӵ�ͼƬ�±�
      int nextIndex = matching_index[sourceIndex][i];

      
      if (adjacent[sourceIndex][nextIndex]) {
        adjacent[sourceIndex][nextIndex] = adjacent[nextIndex][sourceIndex] = false;
        unstitched_idx.push(nextIndex);


        cout << "get Features.\n";
        // kd���������
        vector<POINT_PAIR> src_to_dst_pairs = getPointPairsFromFeatures(features[sourceIndex], features[nextIndex]);
        vector<POINT_PAIR> dst_to_src_pairs = getPointPairsFromFeatures(features[nextIndex], features[sourceIndex]);
        // �����ƥ�䷽��
        if (src_to_dst_pairs.size() > dst_to_src_pairs.size())
          ReplacePairs(src_to_dst_pairs, dst_to_src_pairs);
        else
          ReplacePairs(dst_to_src_pairs, src_to_dst_pairs);

        // RANSAC�㷨
        cout << "RANSAC" << endl;
        vector<int> indices = RANSAC(dst_to_src_pairs);

        // �������ģ�ͣ�����ƽ��λ��
        vector<int> offset = getAvgOffset(src_to_dst_pairs, indices);

        cout << "offset_x " << offset[0] << " offset_y" << offset[1] << endl;

        // ʹ��ƽ��λ����Ϣƴ��ͼƬ
        cur_stitched_img = blending(cur_stitched_img, src_imgs[nextIndex], offset[0], offset[1], offset[2], offset[3]);
        cur_stitched_img.display("mid-process", false);
      }
    }
  }
  return cur_stitched_img;
}