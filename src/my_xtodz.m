function varargout=my_xtodz(x,L,varargin)
% 函数名称：my_xtodz
% 函数功能：载荷挠曲模型
% 输入：L      :标称位置     
%       u      :机翼变形角
% 输出：Rf:子惯位移
dz=0;
if 3==nargin
    u=varargin{1};
    dz=L*sin(u)*cos(u);
    xr=3*(L-L*sin(u)*sin(u))^4;
    qD=dz/xr;
end
if 4==nargin
    u=varargin{1};
    qD=varargin{2};
end
dz=qD*(x.^4-4*L*x.^3+6*L^2*x.^2);
if 1==nargout
    varargout{1}=dz;
end
if 2==nargout
    varargout{1}=dz;
    varargout{2}=qD;
end