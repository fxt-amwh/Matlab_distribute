%% 相控阵天线方向图
%考虑到方便打包，不设置额外封装函数
clc
clear
close all
sita=linspace(-pi/2,pi/2,180);
phi=linspace(-pi,pi,180);
sita0=0;%由控制实现的波束指向
phi0=0;%由控制实现的波束指向
M=33;%x方向阵元数量
N=21;%y方向阵元数量
MList=0:M-1;
NList=0:N-1;
[MList,NList]=meshgrid(MList,NList);
cm=0.01;
lambda=5*cm;%波长
dx=0.6*lambda;%x方向阵元距离
dy=0.6*lambda;%y方向阵元距离
amn=1;%参考单元电流幅度
k=2*pi/lambda;%位置->相位常数
% f=sin((cos(sita).*sin(phi)-1)*(N/2)*pi)./(sin((cos(sita).*sin(phi)-1)*pi/2)*N);
I0=1;%阵子中心电流幅度
n=120*pi;%波阻抗
% f=I0*n/(2*pi)*(sin(pi/2.*cos(sita)))./cos(sita);
f=zeros(length(phi),length(sita));
for ik=1:length(phi)
    for i=1:length(sita)
        f(ik,i)=I0*n/(2*pi)*(sin(pi/2.*cos(sita(i))))./cos(sita(i))*...
        abs(1*sum(sum(amn*...
            exp(1j*k*(...
            MList*dx*(sin(sita(i))*cos(phi(ik))-sin(sita0)*cos(phi0))...
            +NList*dy*(sin(sita(i))*sin(phi(ik))-sin(sita0)*sin(phi0))...
            )))));
    end
end
x=(f.*(ones(length(phi),1)*sin(sita))).*(cos(phi)'*ones(1,length(sita)));
y=(f.*(ones(length(phi),1)*sin(sita))).*(sin(phi)'*ones(1,length(sita)));
z=f.*(ones(length(phi),1)*cos(sita));
figure
surf(x,y,z);
axis equal
title('11*11方阵相控阵方向三维图','Fontsize',15);
xlabel('x','Fontsize',15);ylabel('y','Fontsize',15);zlabel('z','Fontsize',15);
temp=zlim;
zlim([0,temp(2)]);
figure
polar(sita,f(1,:));
title('10*10方阵相控阵方向图');
figure
plot(sita,f(1,:),'b');
hold on;
grid on;
xlabel('theta/radian');
ylabel('amplitude');
title('归一化方向图phi=0');