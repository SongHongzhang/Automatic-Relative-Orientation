%2021/5/26
%本程序实现 Harris+Surf特征点自动提取、匹配，并基于匹配点进行相对定向
clc,clear
tic;
%读入图像，并转换为灰度图
str1='C:\Users\Administrator\Desktop\Photogrammetry\Data\01156.jpg';
str2='C:\Users\Administrator\Desktop\Photogrammetry\Data\01157.jpg';
i1= imread(str1);
I1=rgb2gray(i1);
i2= imread(str2);
I2=rgb2gray(i2);

%harris+surf特征点提取
points1_harr = detectHarrisFeatures(I1,'MinQuality',0.05);
points2_harr = detectHarrisFeatures(I2,'MinQuality',0.05);
points1_surf = detectSURFFeatures(I1,'MetricThreshold',3000);
points2_surf = detectSURFFeatures(I2,'MetricThreshold',3000);

%存储匹配点坐标
resultpairs1 = vertcat(points1_harr.Location,points1_surf.Location);
resultpairs2 = vertcat(points2_harr.Location,points2_surf.Location);
points1=cornerPoints(resultpairs1);
points2=cornerPoints(resultpairs2);

%计算描述向量
[f1, vpts1] = extractFeatures(I1, points1);
[f2, vpts2] = extractFeatures(I2, points2);

%进行匹配
indexPairs = matchFeatures(f1, f2,'Method','NearestNeighborRatio','MaxRatio',0.75) ;
matched_pts1 = vpts1(indexPairs(:, 1));
matched_pts2 = vpts2(indexPairs(:, 2));

%存储匹配点坐标
resultpairs1 = matched_pts1.Location;
resultpairs2 = matched_pts2.Location;

before=size(resultpairs1,1);
%删除数组中重复
result_surf1=unique(resultpairs1,'rows');
result_surf2=unique(resultpairs2,'rows');


%剔除误匹配 (RANSAC)
[hs_ransac, inliers] = estimateFundamentalMatrix(resultpairs1,resultpairs2, 'Method', 'RANSAC','DistanceThreshold',0.1);

XY1=resultpairs1(inliers,:);
XY2=resultpairs2(inliers,:);
%XY1=cell2mat(struct2cell(load('XY1.mat')));
%XY2=cell2mat(struct2cell(load('XY2.mat')));
after=size(XY1,1);

%连续法相对定向
%主距，像平面坐标原点
f=152.720000;
x_01=2528.957;
y_01=2563.180;
x_02=2558.568;
y_02=2559.847;

%同名点像平面坐标
x_1=XY1(:,1)'-x_01;
x_2=XY2(:,1)'-x_02;
y_1=XY1(:,2)'-y_01;
y_2=XY2(:,2)'-y_02;

%相对定向元素初值
phi=0;
omega=0;
kappa=0;
mu=0;
nu=0;

%相对定向点个数
n=length(XY1(:,1));

%改正数初值
d_phi=0;
d_omega=0;
d_kappa=0;
d_mu=0;
d_nu=0;

%限差
limit=1e-6;

%像空间坐标
X1=x_1;
Y1=y_1;
round=0;
for i=1:n
    Z1(i)=-f;
end

%解算法方程
while(1)
    a1=cos(phi)*cos(kappa)-sin(phi)*sin(omega)*sin(kappa);
    a2=-cos(phi)*sin(kappa)-sin(phi)*sin(omega)*cos(kappa);
    a3=-sin(phi)*cos(omega);
    b1=cos(omega)*sin(kappa);
    b2=cos(omega)*cos(kappa);
    b3=-sin(omega);
    c1=sin(phi)*cos(kappa)+cos(phi)*sin(omega)*sin(kappa);
    c2=-sin(phi)*sin(kappa)+cos(phi)*sin(omega)*cos(kappa);
    c3=cos(phi)*cos(omega);
    R=[a1,a2,a3;
        b1,b2,b3;
        c1,c2,c3];
    
    %像空间辅助坐标系
    for i=1:n
        a=(R*[x_2(i);y_2(i);-f])';
        X2(i)=a(1);Y2(i)=a(2);Z2(i)=a(3);
    end
    
    %假定摄影基线
    BX=x_1(1)-x_2(1);
    BY=mu*BX;
    BZ=nu*BX;
    coe= zeros(n,5);
    for i = 1:n
        N =(BX*Z2(i)-BZ*X2(i))/(X1(i)*Z2(i)-Z1(i)*X2(i));
        N_1=(BX*Z1(i)-BZ*X1(i))/(X1(i)*Z2(i)-Z1(i)*X2(i));
        Q=N*Y1(i)-N_1*Y2(i)-BY;
        L(i)=Q;
        coe(i,:)=[-X2(i)*Y2(i)*N_1/Z2(i),-(Z2(i)+Y2(i)^2/Z2(i))*N_1,X2(i)*N_1,BX,-Y2(i)/Z2(i)*BX];
    end
    %权均为1，对角阵省略
    dX=inv(coe'*coe)*coe'*L';
    
    d_phi=dX(1);
    d_omega=dX(2);
    d_kappa=dX(3);
    d_mu=dX(4);
    d_nu=dX(5);
    
    phi=phi+d_phi;
    omega=omega+d_omega;
    kappa=kappa+d_kappa;
    mu=mu+d_mu;
    nu=nu+d_nu;
    round = round+1;
    if(abs(d_phi)<limit && abs(d_omega)<limit && abs(d_kappa)<limit && abs(d_mu)<limit && abs(d_nu)<limit)
        V=coe*dX-L';
        sigma_0=sqrt((V'*V)/(n-5))
        break;
    end
end
Q_xx=inv(coe'*coe);
m_mu=Q_xx(1,1);
m_nu=Q_xx(2,2);
m_phi=Q_xx(3,3);
m_omega=Q_xx(4,4);
m_kappa=Q_xx(5,5);
end

% 显示有误匹配的情况
figure
showMatchedFeatures(i1,i2,double(resultpairs1),double(resultpairs2),'montage');
str=sprintf('Matched points');
title(str,'fontname','Times New Roman','FontSize',12);
legend('matched points 1','matched points 2');


% 显示没有误匹配的情况
figure
showMatchedFeatures(i1,i2,resultpairs1(inliers,:),resultpairs2(inliers,:),'montage');
str=sprintf('Matched inlier points');
title(str,'fontname','Times New Roman','FontSize',12);
legend('matched points 1','matched points 2');

fid=fopen('Result.txt','w');%写入文件路径
fprintf(fid , '基于 Harris 与 SURF 特征点的影像自动相对定向结果\n\n');
fprintf(fid,'特征匹配\n');
fprintf(fid,'粗匹配得到的同名特征点对数：');
fprintf(fid,'%d\n',before);
fprintf(fid,'RANSAC 去除误匹配后得到的同名特征点对数：');
fprintf(fid,'%d\n\n',after);
fprintf(fid,'相对定向\n');
fprintf(fid,'参与相对定向的同名特征点对数：');
fprintf(fid,'%d\n',n);
fprintf(fid,'相对定向参数\n');
fprintf(fid,'μ=');
fprintf(fid,'%f',mu);
fprintf(fid,'±');
fprintf(fid,'%.16f\n',m_mu);
fprintf(fid,'ν=');
fprintf(fid,'%f',nu);
fprintf(fid,'±');
fprintf(fid,'%.16f\n',m_nu);
fprintf(fid,'φ=');
fprintf(fid,'%f',phi);
fprintf(fid,'±');
fprintf(fid,'%.16f\n',m_phi);
fprintf(fid,'ω=');
fprintf(fid,'%f',omega);
fprintf(fid,'±');
fprintf(fid,'%.16f\n',m_omega);
fprintf(fid,'κ=');
fprintf(fid,'%f',kappa);
fprintf(fid,'±');
fprintf(fid,'%.16f\n',m_kappa);
fprintf(fid,'单位权中误差：σ0=');
fprintf(fid,'%f\n',sigma_0);
fclose(fid)
toc