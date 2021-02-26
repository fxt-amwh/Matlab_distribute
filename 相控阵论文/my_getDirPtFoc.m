function ROpf=my_getDirPtFoc(RO_base,varargin)
% 函数名称：my_getDirPtFoc
% 函数功能：由子惯输出进行相对导航滤波
% 输入：RO_base      :雷达基本结构体         
% 输出：ROpf:雷达方向图
% phi=RO_base.phi;
% sita=RO_base.sita;
% RO_base.lambda;%波长
% RO_base.M;%x方向阵元个数
% RO_base.N;%y方向阵元个数
% RO_base.sita0;%波束控制指向俯仰角
% RO_base.phi0;%波束控制指向方位角
% RO_base.ddsita;%俯仰角误差矩阵
% RO_base.ddphi;%方位角误差矩阵
% RO_base.ddx;%x方向位置误差矩阵
% RO_base.ddy;%y方向位置误差矩阵
% RO_base.ddz;%z方向位置误差矩阵
phi=RO_base.phi;
sita=RO_base.sita;
lambda=RO_base.lambda;%波长
M=RO_base.M;
N=RO_base.N;
ddsita=RO_base.ddsita;
ddphi=RO_base.ddphi;
sita0=RO_base.sita0;%由控制实现的波束指向
phi0=RO_base.phi0;%由控制实现的波束指向
ddx=RO_base.ddx;
ddy=RO_base.ddy;
ddz=RO_base.ddz;


MList=0:M-1;
NList=0:N-1;
[MList,NList]=meshgrid(MList,NList);
% dx=lambda/2;%x方向阵元距离
% dy=lambda/2;%y方向阵元距离
dx=RO_base.dx;%x方向阵元距离
dy=RO_base.dy;%y方向阵元距离
amn=1;%参考单元电流幅度
k=2*pi/lambda;%位置->相位常数
I0=1;%阵子中心电流幅度
n=120*pi;%波阻抗
    
ROpf=zeros(length(phi),length(sita));
for ik=1:length(phi)
    for i=1:length(sita)
        ROpf(ik,i)=I0*n/(2*pi)*...
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