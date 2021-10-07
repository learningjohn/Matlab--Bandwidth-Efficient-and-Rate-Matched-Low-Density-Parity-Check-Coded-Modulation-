# 尝试用matlab实现概率成型算法
Matlab实现论文“Bandwidth Efficient and Rate-Matched Low-Density Parity-Check Coded Modulation”中的概率整形算法

第一次学习使用Matlab中的工程功能，并用git进行维护

10月7号的版本，实现内容:

1)、完成了论文中提出的二分法和黄金分割法计算的最优![图片](https://user-images.githubusercontent.com/48309583/136379048-b13741fb-fa05-4ef0-99eb-e9a6b5363b56.png)值
并计算最优的M-B分布P。
（main函数）

2)、产生满足最优的M-B分布P下的随机QAM符号 与 满足均匀分布产生的QAM符号，这两个在AWGN下进行性能比较。前面P分布的随机QAM符号通过CCDM产生。结果表明最优分布下的QAM符号在容量有优势的部分其误码率性能也有提升。
（test1函数）

3）、比较随机信号和图像信号的信道互信息量，图像信号先转换为比特，在通过比特转换为QAM符号。图像信号的比特提取有两种方式，一种是对各个像素点进行排列，一起从高位取一直取到低位。另一种是将像素点转换为比特，8比特依次排列取。两种取法得到的互信息量不同。高位到低位的取法互信息量更低。
