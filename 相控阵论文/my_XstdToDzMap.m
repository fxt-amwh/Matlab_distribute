function [ddx,ddz]=my_XstdToDzMap(xstd,L,u,varargin)
% 函数名称：my_xtodz
% 函数功能：载荷挠曲模型 标准位置在挠曲后的位置
% 输入：xstd      :求取位置    
%      L      :机翼翼尖位置
%      u      :机翼变形角
% 输出：varargout{1}:z轴位移
%     ：varargout{2}实际位置误差
err=1e-5;
flag=0;
if 4 == nargin
    err=varargin{1};
elseif 5 == nargin
    err=varargin{1};
    flag=varargin{2};
end
[N,M]=size(xstd);
ddx=zeros(N,M);
ddz=zeros(N,M);
for n=1:N
    for m=1:M
        [dz,dx]=my_XstdToDz(xstd(n,m),L,u,err,flag);
        ddx(n,m)=dx;ddz(n,m)=dz;
    end
end