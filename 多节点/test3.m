%% 载荷挠曲模型测试
% 版本时间：2020/12/20 
% 设计初衷：载荷挠曲模型测试
close all;
x=linspace(0,5,100);
xstd=linspace(0,5,100);
L=5;
u=5*pi/180;
dz=my_xtodz(x,L,u);
[x1,dz1]=my_xtodz1(xstd,L,u);
figure
plot(x,dz);
hold on;
plot(x1,dz1);
title("载荷挠曲模型")
ylim([0,0.5]);
figure
plot(x,dz-dz1);
title("载荷挠曲误差模型")