%% ��������߷���ͼ
%���ǵ��������������ö����װ���� 14�������
%�汾ʱ�䣺2020/01/21
%��Ƴ��ԣ�����ͼfor����
clc
clear
close all

global lambda
global M
global N
sita=linspace(-pi/2,pi/2,180);
phi=linspace(0,pi,180);
cm=0.01;
M=33;%x������Ԫ����
N=21;%y������Ԫ���� ��õ���M
lambda=5*cm;%����
dx=lambda*0.6;%x������Ԫ����
dy=lambda*0.6;%y������Ԫ����
xm=floor(M/2)*dx;%x����������
ym=floor(N/2)*dy;%y����������
% load NSINSFR3;
% % umax=NSINSFR.attnewall_true(2,1:2:7000)/2;%��λ������������ʵֵ �غ�ģ��
% % umax_measure=NSINSFR.attnewall(2,1:2:7000)/2;%��λ��������������ֵ �غ�ģ��
% umax=NSINSFR.attnewall_true(2,125)/2;%��λ������������ʵֵ �غ�ģ��
% umax_measure=NSINSFR.attnewall(2,125)/2;%��λ��������������ֵ �غ�ģ��
% ddx=0*0.001*randn(N,M);
% ddy=0*0.001*randn(N,M);
% ddz=0*0.001*randn(N,M);
% xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%��Ƥ���İ�װ�� 2m ��
% [d_x,dz]=my_fixdefload(xdir+ddx,1,umax(loopi),1);%�����غ�����ģ�ͱ�׼�� ��׼ģ��1
% [dx_measure,dz_measure]=my_fixdefload(xdir+ddx,1,umax_measure(loopi),1);%�����غ�����ģ�Ͳ�����
f=my_getDirPt(sita,phi,0,0,0);%��λ������ͼ����
x=(f.*(ones(length(phi),1)*sin(sita))).*(cos(phi)'*ones(1,length(sita)));
y=(f.*(ones(length(phi),1)*sin(sita))).*(sin(phi)'*ones(1,length(sita)));
z=f.*(ones(length(phi),1)*cos(sita));
figure
surf(x,y,z);
axis equal
title(sprintf('%d*%d�������������άͼ',N,M),'Fontsize',15);
xlabel('x','Fontsize',15);ylabel('y','Fontsize',15);zlabel('z','Fontsize',15);
temp=zlim;
zlim([0,temp(2)]);
figure
polar(sita,f(1,:).*sin(phi)/max(f(1,:)));
title(sprintf('%d*%d���������H�淽��ͼ,d=2*lamda,��λ=�ͺ�',N,M));
figure
plot(sita,f(1,:)/max(f(1,:)),'b');
hold on;
grid on;
xlabel('theta/radian');
ylabel('amplitude');
title('��һ������ͼphi=0');