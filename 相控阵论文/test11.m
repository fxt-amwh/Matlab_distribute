%% ��������߷���ͼ
%���ǵ��������������ö����װ����
clc
clear
close all
sita=linspace(-pi/2,pi/2,180);
phi=linspace(-pi,pi,180);
sita0=0;%�ɿ���ʵ�ֵĲ���ָ��
phi0=0;%�ɿ���ʵ�ֵĲ���ָ��
M=33;%x������Ԫ����
N=21;%y������Ԫ����
MList=0:M-1;
NList=0:N-1;
[MList,NList]=meshgrid(MList,NList);
cm=0.01;
lambda=5*cm;%����
dx=0.6*lambda;%x������Ԫ����
dy=0.6*lambda;%y������Ԫ����
amn=1;%�ο���Ԫ��������
k=2*pi/lambda;%λ��->��λ����
% f=sin((cos(sita).*sin(phi)-1)*(N/2)*pi)./(sin((cos(sita).*sin(phi)-1)*pi/2)*N);
I0=1;%�������ĵ�������
n=120*pi;%���迹
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
title('11*11�������������άͼ','Fontsize',15);
xlabel('x','Fontsize',15);ylabel('y','Fontsize',15);zlabel('z','Fontsize',15);
temp=zlim;
zlim([0,temp(2)]);
figure
polar(sita,f(1,:));
title('10*10�����������ͼ');
figure
plot(sita,f(1,:),'b');
hold on;
grid on;
xlabel('theta/radian');
ylabel('amplitude');
title('��һ������ͼphi=0');