## 感谢评价本次作业！



##  文件结构

1. data1是用于测试均衡化的数据集
2. data2是用于测试color transform的数据集
3. images中存放实验截图



## 运行方式

```bash
$ g++ Equalization.cpp -o test -O2 -lgdi32 -std=c++11
$ test 图片1 图片2 op

(op为0是均衡化，第二个参数无效，随便输就可以；op为1是color transform，图片1是原图，图片2是目标图片)
```



### 如有任何问题，请与我联系：jacky14.liang@gmail.com，非常感谢！