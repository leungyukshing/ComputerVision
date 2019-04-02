from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score
from sklearn.externals import joblib
import tensorflow.examples.tutorials.mnist.input_data as input_data
import time
import os
from skimage import io,data,transform
from datetime import datetime
import numpy as np

MNIST_SIZE = 28

def translate(image_path):
   #读入图片并变成灰色
  img = io.imread(image_path, as_grey=True)
  #缩小到28*28
  translated_img = transform.resize(img, (MNIST_SIZE, MNIST_SIZE))
  #变成1*784的一维数组
  flatten_img = np.reshape(translated_img, 784)
  #print(flatten_img)
  #mnist数据集中1代表黑，0代表白
  #result = np.array(255* [flatten_img])
  result = flatten_img * 255
  #print(result)
  #返回该图的所代表的向量
  return result

def train_model(clf_rf, save_dir):
  StartTime = time.clock()
  # Training
  clf_rf.fit(batch_x,batch_y)

  # Save Model
  joblib.dump(clf_rf, save_dir)

  EndTime = time.clock()
  print('Total time %.2f s' % (EndTime - StartTime))


if __name__ == '__main__':
  print("Start")
  for index in range(0, 18):
    mnist = translate("./result/SingleNumImg/16340311.bmp/2_" + str(index) + ".bmp");
    batch_size = 50000
     
    #print("start AdaBoosting")
    save_dir = './model/clf.m'
    i = 10

    if(os.path.isfile(save_dir)):
      clf_rf = joblib.load(save_dir)
    else:
      clf_rf = AdaBoostClassifier(DecisionTreeClassifier(max_depth=5, min_samples_split=5, min_samples_leaf=5), n_estimators=i, learning_rate=0.05, algorithm='SAMME.R')
      train_model(clf_rf, save_dir)

    y_pred_rf = clf_rf.predict([mnist])
    print("Result: %d", y_pred_rf)
