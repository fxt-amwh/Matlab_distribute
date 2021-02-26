clc
clear
close all
sita=linspace(0,2*pi,100);
phi=linspace(0,pi,100);
N=14;
f=sin((cos(sita).*sin(phi)-1)*(N/2)*pi)./(sin((cos(sita).*sin(phi)-1)*pi/2)*N);
% I0=1;%阵子中心电流幅度
% n=120*pi;%波阻抗
% f=I0*n/(2*pi)*(sin(pi/2.*cos(sita)))./cos(sita);
figure
polar(sita,f.*sin(phi));
title('14元端射振的H面方向图,d=lamda/2,相位=滞后');
x=(f.*sin(sita))'*cos(phi);
y=(f.*sin(sita))'*sin(phi);
z=(f.*cos(sita))'*ones(size(phi));
figure
surf(x,y,z);
axis equal
title('14元阵列方向三维图')
xlabel('x');ylabel('y');zlabel('z');