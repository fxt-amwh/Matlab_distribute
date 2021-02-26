function ROpf=my_getDirPtFoc(RO_base,varargin)
% �������ƣ�my_getDirPtFoc
% �������ܣ����ӹ����������Ե����˲�
% ���룺RO_base      :�״�����ṹ��         
% �����ROpf:�״﷽��ͼ
% phi=RO_base.phi;
% sita=RO_base.sita;
% RO_base.lambda;%����
% RO_base.M;%x������Ԫ����
% RO_base.N;%y������Ԫ����
% RO_base.sita0;%��������ָ������
% RO_base.phi0;%��������ָ��λ��
% RO_base.ddsita;%������������
% RO_base.ddphi;%��λ��������
% RO_base.ddx;%x����λ��������
% RO_base.ddy;%y����λ��������
% RO_base.ddz;%z����λ��������
phi=RO_base.phi;
sita=RO_base.sita;
lambda=RO_base.lambda;%����
M=RO_base.M;
N=RO_base.N;
ddsita=RO_base.ddsita;
ddphi=RO_base.ddphi;
sita0=RO_base.sita0;%�ɿ���ʵ�ֵĲ���ָ��
phi0=RO_base.phi0;%�ɿ���ʵ�ֵĲ���ָ��
ddx=RO_base.ddx;
ddy=RO_base.ddy;
ddz=RO_base.ddz;


MList=0:M-1;
NList=0:N-1;
[MList,NList]=meshgrid(MList,NList);
% dx=lambda/2;%x������Ԫ����
% dy=lambda/2;%y������Ԫ����
dx=RO_base.dx;%x������Ԫ����
dy=RO_base.dy;%y������Ԫ����
amn=1;%�ο���Ԫ��������
k=2*pi/lambda;%λ��->��λ����
I0=1;%�������ĵ�������
n=120*pi;%���迹
    
ROpf=zeros(length(phi),length(sita));
for ik=1:length(phi)
    for i=1:length(sita)
        ROpf(ik,i)=I0*n/(2*pi)*...
        abs(1*sum(sum(amn*(sin(pi/2.*cos(sita(i)+ddsita)))./cos(sita(i)+ddsita).*...
            exp(1j*k*(...
            MList.*dx.*(sin(sita(i)+ddsita).*cos(phi(ik)+ddphi)-sin(sita0)*cos(phi0))...
            +NList.*dy.*(sin(sita(i)+ddsita).*sin(phi(ik)+ddphi)-sin(sita0)*sin(phi0))...
            +(ddx-ddx(ceil(size(ddx,1)/2),ceil(size(ddx,2)/2))).*(sin(sita(i)+ddsita).*cos(phi(ik)+ddphi)-sin(sita0)*cos(phi0))...%��������Ϊ���
            +(ddy-ddy(ceil(size(ddx,1)/2),ceil(size(ddx,2)/2))).*(sin(sita(i)+ddsita).*sin(phi(ik)+ddphi)-sin(sita0)*sin(phi0))...
            +(ddz-ddz(ceil(size(ddx,1)/2),ceil(size(ddx,2)/2))).*cos(sita(i)+ddsita)...
            )))));
    end
end