## 文件说明

1. `canny.cpp`和`canny.h`是我实现封装的Canny算法。
2. `canny_source.c`是老师提供的Canny算法代码。
3. `CImg.h`是实现过程中用到的CImg库。
4. 其余`.png`文件是四个不同测试数据的测试结果截图。

## 运行方式

在本目录下打开命令行，输入：

```bash
g++ canny.cpp -o test -std=c++11 -O2 -lgdi32
test filename sigma tlow thigh
```



### 感谢您对我的作业进行评价！

