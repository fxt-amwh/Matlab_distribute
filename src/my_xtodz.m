function varargout=my_xtodz(x,L,varargin)
% �������ƣ�my_xtodz
% �������ܣ��غ�����ģ��
% ���룺L      :���λ��     
%       u      :������ν�
% �����Rf:�ӹ�λ��
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