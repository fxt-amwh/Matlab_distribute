function varargout=my_XstdToDz(xstd,L,u,varargin)
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
if flag==2
    qD=u;
else
    [~,qD]=my_xtodz(0,L,u);%根据翼尖位置、挠曲角求取载荷系数
end
fw=@(x)sqrt(1+qD*(4*x.^3-12*L*x.^2+12*L^2*x));
xloop=xstd;
dx=1;
i=0;
while abs(dx) > err
    dx=integral(fw,0,xloop)-xstd;
    xloop=xloop-dx;
    i=i+1;
end
if 1==nargout
    varargout{1}=my_xtodz(xloop,L,u,qD);
end
if 2==nargout
    varargout{1}=my_xtodz(xloop,L,u,qD);
    varargout{2}=xloop-xstd;
end
if 3==nargout
    varargout{1}=my_xtodz(xloop,L,u,qD);
    varargout{2}=xloop-xstd;
    varargout{3}=i;
end