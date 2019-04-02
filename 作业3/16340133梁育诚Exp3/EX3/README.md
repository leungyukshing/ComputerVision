## 运行方式

编译：在命令行输入：

```bash
g++ hough.cpp -o test -O2 -lgdi32 -std=c++1
```

运行：在命令行输入：

```bash
test ./数据集/文件名 sigma tlow thigh minRadius maxRadius
```



## 文件结构说明

`set1`和`set2`是存放测试数据的文件夹。