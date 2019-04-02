## 运行方式

### 第一部分

编译：在命令行输入：

```bash
g++ main.cpp Segmentation -o test -O2 -lgdi32 -std=c++1
```

运行：在命令行输入：

```java
test
```

---

### 第二部分

训练模型

```bash
python adaboost.py
```

若要更改estimator数量，进入`adaboost.py`中修改`i`的值即可。

测试数据

```bash
python test.py
```

