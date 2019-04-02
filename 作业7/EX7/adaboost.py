from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score
from sklearn.externals import joblib
import tensorflow.examples.tutorials.mnist.input_data as input_data
import time
import os
from datetime import datetime

def train_model(clf_rf, save_dir):
  StartTime = time.clock()
  # Training
  clf_rf.fit(batch_x,batch_y)

  EndTime = time.clock()
  
  # Save Model
  joblib.dump(clf_rf, save_dir)

  print('Total time %.2f s' % (EndTime - StartTime))


if __name__ == '__main__':
  data_dir = './MNIST_data/'
  mnist = input_data.read_data_sets(data_dir,one_hot=False)
  batch_size = 50000
  batch_x,batch_y = mnist.train.next_batch(batch_size)
  test_x = mnist.test.images[:10000]
  test_y = mnist.test.labels[:10000]
   
  print("start AdaBoosting")
  i = 30
  save_dir = './model/clf_rf_' + str(i) +'.m' 
  

  if(os.path.isfile(save_dir)):
    clf_rf = joblib.load(save_dir)
  else:
    clf_rf = AdaBoostClassifier(DecisionTreeClassifier(max_depth=5, min_samples_split=5, min_samples_leaf=5), n_estimators=i, learning_rate=0.05, algorithm='SAMME.R')
    train_model(clf_rf, save_dir)

  y_pred_rf = clf_rf.predict(test_x)
  print("Result: %d", y_pred_rf)
  acc_rf = accuracy_score(test_y,y_pred_rf)
  print("%s n_estimators = %d, random forest accuracy:%f" % (datetime.now(), i, acc_rf))
