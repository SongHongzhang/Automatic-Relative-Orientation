# Automatic-Relative-Orientation
&emsp;&emsp;数字摄影测量实习之自动相对定向. 先提取 Harris 特征与 SURF 特征进行特征匹配，再使用 RANSAC 算法去除粗差得到同名特征点，最后利用这些同名点进行相对定向的解算. 

- 使用 `tic` 和 `toc` 可以得到程序运行的时间，约 18 秒
- 在相对定向的过程中没有对视差较大的点进行逐点剔除
- 得到同名特征点的过程使用内置的函数完成
- 并没有进行相对定向精度的控制，偶尔可以达到 1/3 像素的精度，最佳（误差达标）的一组数据点见 *XY1.mat*，*XY2.mat*，使用 `XY1=cell2mat(struct2cell(load('XY1.mat')));` 导入这些数据

&emsp;&emsp;加油
