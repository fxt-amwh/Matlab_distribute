function f=my_getDirPt(sita,phi,ddx,ddy,ddz,varargin)
global lambda;%����
global M;
global N;
ddsita=0;
ddphi=ddsita;
if nargin>5
    if nargin==6
        ddsita=varargin{1};
    elseif nargin==7
        ddsita=varargin{1};
        ddphi=varargin{2};
    else
        
    end
else
    
end
sita0=0;%�ɿ���ʵ�ֵĲ���ָ��
phi0=0;%�ɿ���ʵ�ֵĲ���ָ��

MList=0:M-1;
NList=0:N-1;
[MList,NList]=meshgrid(MList,NList);
dx=lambda*0.6;%x������Ԫ����
dy=lambda*0.6;%y������Ԫ����
amn=1;%�ο���Ԫ��������
k=2*pi/lambda;%λ��->��λ����
I0=1;%�������ĵ�������
n=120*pi;%���迹
    
f=zeros(length(phi),length(sita));
for ik=1:length(phi)
    for i=1:length(sita)
        f(ik,i)=I0*n/(2*pi)*...
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

