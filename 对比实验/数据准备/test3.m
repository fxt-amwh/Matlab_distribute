%% �غ�����ģ�Ͳ���
% �汾ʱ�䣺2020/12/20 
% ��Ƴ��ԣ��غ�����ģ�Ͳ���
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
title("�غ�����ģ��")
ylim([0,0.5]);
figure
plot(x,dz-dz1);
title("�غ��������ģ��")