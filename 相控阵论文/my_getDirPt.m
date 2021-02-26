function f=my_getDirPt(sita,phi,ddx,ddy,ddz,varargin)
global lambda;%波长
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
sita0=0;%由控制实现的波束指向
phi0=0;%由控制实现的波束指向

MList=0:M-1;
NList=0:N-1;
[MList,NList]=meshgrid(MList,NList);
dx=lambda*0.6;%x方向阵元距离
dy=lambda*0.6;%y方向阵元距离
amn=1;%参考单元电流幅度
k=2*pi/lambda;%位置->相位常数
I0=1;%阵子中心电流幅度
n=120*pi;%波阻抗
    
f=zeros(length(phi),length(sita));
for ik=1:length(phi)
    for i=1:length(sita)
        f(ik,i)=I0*n/(2*pi)*...
        abs(1*sum(sum(amn*(sin(pi/2.*cos(sita(i)+ddsita)))./cos(sita(i)+ddsita).*...
            exp(1j*k*(...
            MList.*dx.*(sin(sita(i)+ddsita).*cos(phi(ik)+ddphi)-sin(sita0)*cos(phi0))...
            +NList.*dy.*(sin(sita(i)+ddsita).*sin(phi(ik)+ddphi)-sin(sita0)*sin(phi0))...
            +(ddx-ddx(ceil(size(ddx,1)/2),ceil(size(ddx,2)/2))).*(sin(sita(i)+ddsita).*cos(phi(ik)+ddphi)-sin(sita0)*cos(phi0))...%以下三项为误差
            +(ddy-ddy(ceil(size(ddx,1)/2),ceil(size(ddx,2)/2))).*(sin(sita(i)+ddsita).*sin(phi(ik)+ddphi)-sin(sita0)*sin(phi0))...
            +(ddz-ddz(ceil(size(ddx,1)/2),ceil(size(ddx,2)/2))).*cos(sita(i)+ddsita)...
            )))));
    end
end

