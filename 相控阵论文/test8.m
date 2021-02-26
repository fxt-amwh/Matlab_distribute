%% ��������߷���ͼ
%���ǵ��������������ö����װ����
%�汾ʱ�䣺2020/1/23
%��Ƴ��ԣ�����ͼfor���� �Ա�
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
M=33;%x������Ԫ����
N=21;%y������Ԫ���� ��õ���M
lambda=5*cm;%����
dx=lambda*0.6;%x������Ԫ����
dy=lambda*0.6;%y������Ԫ����
xm=floor(M/2)*dx;%x����������
ym=floor(N/2)*dy;%y����������
load NSINSFR3;
% umax=NSINSFR.attnewall_true(2,1:2:7000)/2;%��λ������������ʵֵ �غ�ģ��
% umax_measure=NSINSFR.attnewall(2,1:2:7000)/2;%��λ��������������ֵ �غ�ģ��
umax=NSINSFR.attnewall_true(2,125)/2;%��λ������������ʵֵ �غ�ģ��
umax_measure=NSINSFR.attnewall(2,125)/2;%��λ��������������ֵ �غ�ģ��
ddx=0.001*randn(N,M);
ddy=0.001*randn(N,M);
ddz=0.001*randn(N,M);
% ���µ������
xdir=0.5 + (-floor(M/2):1:floor(M/2))*dx;%��Ƥ���İ�װ�� 0.5m ��
[d_x,dz]=my_fixdefload(xdir+ddx,2,10*umax,1);%�����غ�����ģ�ͱ�׼�� ��׼ģ��1
[dx_measure,dz_measure]=my_fixdefload(xdir+ddx,2,10*umax_measure,1);%�����غ�����ģ�Ͳ�����

f_now=my_getDirPt(sita,phi,d_x-dx_measure,ddy,dz+ddz-dz_measure);%��������ͼ����
f_5=my_getDirPt(sita,phi,d_x,ddy,dz+ddz);%5������ָ������ͼ����
f=my_getDirPt(sita,0,0,0,0);%��λ�����
% x=(f_now.*(ones(length(phi),1)*sin(sita))).*(cos(phi)'*ones(1,length(sita)));
% y=(f_now.*(ones(length(phi),1)*sin(sita))).*(sin(phi)'*ones(1,length(sita)));
% z=f_now.*(ones(length(phi),1)*cos(sita));
% figure
% surf(x,y,z);
% axis equal
% title('���������άͼ','Fontsize',15);
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
xlabel('theta/��');
ylabel('amplitude');
ylim([0,1])
xlim([-90,90])
title('��һ������ͼphi=0��,����ָ��10��');
legend("����ǰ","������",'Fontsize',15);
%%
figure
% load f_5;
polarplot(sita,f_5(1,:).*sin(phi)/max(f),'LineWidth',1);
hold on;
polarplot(sita,f_now(1,:).*sin(phi)/max(f),'LineWidth',1);
title('�����H�淽��ͼ,d=2*lamda,��λ=�ͺ�,����ָ��10��','Fontsize',15);
rlim([0,1])
legend("����ǰ","������",'Fontsize',15);