clc
clear
close all
sita=linspace(0,2*pi,100);
phi=linspace(0,pi,100);
N=14;
f=sin((cos(sita).*sin(phi)-1)*(N/2)*pi)./(sin((cos(sita).*sin(phi)-1)*pi/2)*N);
% I0=1;%�������ĵ�������
% n=120*pi;%���迹
% f=I0*n/(2*pi)*(sin(pi/2.*cos(sita)))./cos(sita);
figure
polar(sita,f.*sin(phi));
title('14Ԫ�������H�淽��ͼ,d=lamda/2,��λ=�ͺ�');
x=(f.*sin(sita))'*cos(phi);
y=(f.*sin(sita))'*sin(phi);
z=(f.*cos(sita))'*ones(size(phi));
figure
surf(x,y,z);
axis equal
title('14Ԫ���з�����άͼ')
xlabel('x');ylabel('y');zlabel('z');