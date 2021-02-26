%% 相控阵天线方向图
%考虑到方便打包，不设置额外封装函数
%版本时间：2020/1/23
%设计初衷：方向图for论文 对比
clc
clear
close all
global lambda
global M
global N
sita=linspace(-pi/2,pi/2,1800);
phi=linspace(0,pi,1800);
% phi=pi/2;
cm=0.01;
M=33;%x方向阵元数量
N=21;%y方向阵元数量 最好等于M
lambda=5*cm;%波长
dx=lambda*0.6;%x方向阵元距离
dy=lambda*0.6;%y方向阵元距离
xm=floor(M/2)*dx;%x方向最大距离
ym=floor(N/2)*dy;%y方向最大距离
load NSINSFR3;
% umax=NSINSFR.attnewall_true(2,1:2:7000)/2;%单位机翼处挠曲角真实值 载荷模型
% umax_measure=NSINSFR.attnewall(2,1:2:7000)/2;%单位机翼处挠曲角量测值 载荷模型
umax=NSINSFR.attnewall_true(2,125)/2;%单位机翼处挠曲角真实值 载荷模型
umax_measure=NSINSFR.attnewall(2,125)/2;%单位机翼处挠曲角量测值 载荷模型
ddx=0.001*randn(N,M);
ddy=0.001*randn(N,M);
ddz=0.001*randn(N,M);
% 最新的相控阵
xdir=0.5 + (-floor(M/2):1:floor(M/2))*dx;%蒙皮中心安装在 0.5m 处
[d_x,dz]=my_fixdefload(xdir+ddx,2,10*umax,1);%机翼载荷挠曲模型标准量 精准模型1
[dx_measure,dz_measure]=my_fixdefload(xdir+ddx,2,10*umax_measure,1);%机翼载荷挠曲模型补偿量

f_now=my_getDirPt(sita,phi,d_x-dx_measure,ddy,dz+ddz-dz_measure);%误差补偿后方向图函数
f_5=my_getDirPt(sita,phi,d_x,ddy,dz+ddz);%5°挠曲指数方向图函数
f=my_getDirPt(sita,0,0,0,0);%无位置误差
% x=(f_now.*(ones(length(phi),1)*sin(sita))).*(cos(phi)'*ones(1,length(sita)));
% y=(f_now.*(ones(length(phi),1)*sin(sita))).*(sin(phi)'*ones(1,length(sita)));
% z=f_now.*(ones(length(phi),1)*cos(sita));
% figure
% surf(x,y,z);
% axis equal
% title('相控阵方向三维图','Fontsize',15);
% xlabel('x','Fontsize',15);ylabel('y','Fontsize',15);zlabel('z','Fontsize',15);
% temp=zlim;
% zlim([0,temp(2)]);
%%
figure
plot(sita*180/pi,f_5(1,:)/max(f));
hold on;
plot(sita*180/pi,f_now(1,:)/max(f));
hold on;
grid on;
xlabel('theta/°');
ylabel('amplitude');
ylim([0,1])
xlim([-90,90])
title('归一化方向图phi=0°,挠曲指数10°');
legend("补偿前","补偿后",'Fontsize',15);
%%
figure
% load f_5;
polarplot(sita,f_5(1,:).*sin(phi)/max(f),'LineWidth',1);
hold on;
polarplot(sita,f_now(1,:).*sin(phi)/max(f),'LineWidth',1);
title('相控阵H面方向图,d=2*lamda,相位=滞后,挠曲指数10°','Fontsize',15);
rlim([0,1])
legend("补偿前","补偿后",'Fontsize',15);