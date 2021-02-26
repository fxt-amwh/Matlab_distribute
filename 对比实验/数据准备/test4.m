%% 载荷挠曲模型测试
% 版本时间：2020/12/20 
% 设计初衷：积分测试
close all;
x=linspace(0,5,100);
xstd=linspace(0,5,100);
L=5;
u=5*pi/180;
[dz,qD]=my_xtodz(x,L,u);
% dz=qD*(x.^4-4*L*x.^3+6*L^2*x.^2);
fw=@(x)sqrt(1+qD*(4*x.^3-12*L*x.^2+12*L^2*x))
ac=@(x)sin(x)./x
s=quad(ac,pi/4,pi/2)
s=integral(fw,0,3)
xstd=1;
[dz,dx,num]=my_XstdToDz(xstd,L,u)
s=integral(fw,0,xstd+dx)