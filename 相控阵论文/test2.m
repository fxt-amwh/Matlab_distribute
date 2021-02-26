%% 相控阵天线方向图
%考虑到方便打包，不设置额外封装函数 14所建议后
%版本时间：2020/01/21
%设计初衷：方向图for论文
clc
clear
close all

global lambda
global M
global N
sita=linspace(-pi/2,pi/2,180);
phi=linspace(0,pi,180);
cm=0.01;
M=33;%x方向阵元数量
N=21;%y方向阵元数量 最好等于M
lambda=5*cm;%波长
dx=lambda*0.6;%x方向阵元距离
dy=lambda*0.6;%y方向阵元距离
xm=floor(M/2)*dx;%x方向最大距离
ym=floor(N/2)*dy;%y方向最大距离
% load NSINSFR3;
% % umax=NSINSFR.attnewall_true(2,1:2:7000)/2;%单位机翼处挠曲角真实值 载荷模型
% % umax_measure=NSINSFR.attnewall(2,1:2:7000)/2;%单位机翼处挠曲角量测值 载荷模型
% umax=NSINSFR.attnewall_true(2,125)/2;%单位机翼处挠曲角真实值 载荷模型
% umax_measure=NSINSFR.attnewall(2,125)/2;%单位机翼处挠曲角量测值 载荷模型
% ddx=0*0.001*randn(N,M);
% ddy=0*0.001*randn(N,M);
% ddz=0*0.001*randn(N,M);
% xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%蒙皮中心安装在 2m 处
% [d_x,dz]=my_fixdefload(xdir+ddx,1,umax(loopi),1);%机翼载荷挠曲模型标准量 精准模型1
% [dx_measure,dz_measure]=my_fixdefload(xdir+ddx,1,umax_measure(loopi),1);%机翼载荷挠曲模型补偿量
f=my_getDirPt(sita,phi,0,0,0);%无位置误差方向图函数
x=(f.*(ones(length(phi),1)*sin(sita))).*(cos(phi)'*ones(1,length(sita)));
y=(f.*(ones(length(phi),1)*sin(sita))).*(sin(phi)'*ones(1,length(sita)));
z=f.*(ones(length(phi),1)*cos(sita));
figure
surf(x,y,z);
axis equal
title(sprintf('%d*%d方阵相控阵方向三维图',N,M),'Fontsize',15);
xlabel('x','Fontsize',15);ylabel('y','Fontsize',15);zlabel('z','Fontsize',15);
temp=zlim;
zlim([0,temp(2)]);
figure
polar(sita,f(1,:).*sin(phi)/max(f(1,:)));
title(sprintf('%d*%d方阵相控阵H面方向图,d=2*lamda,相位=滞后',N,M));
figure
plot(sita,f(1,:)/max(f(1,:)),'b');
hold on;
grid on;
xlabel('theta/radian');
ylabel('amplitude');
title('归一化方向图phi=0');