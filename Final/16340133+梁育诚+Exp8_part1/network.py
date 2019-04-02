from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score
from sklearn.externals import joblib
from sklearn import svm
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
   
  print("start SVM")
  save_dir = './model/clfm' 
  

  if(os.path.isfile(save_dir)):
    clf_rf = joblib.load(save_dir)
  else:
    clf_rf = MLPClassifier(hidden_layer_sizes=(400, 200), activation='logistic', 
        solver='sgd', learning_rate_init=0.001, max_iter=400, verbose = True)
    train_model(clf_rf, save_dir)

  y_pred_rf = clf_rf.predict(test_x)
  print("Result: %d", y_pred_rf)
  acc_rf = accuracy_score(test_y,y_pred_rf)
  print("%s , random forest accuracy:%f" % (datetime.now(), acc_rf))
