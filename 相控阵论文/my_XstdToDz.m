function varargout=my_XstdToDz(xstd,L,u,varargin)
% �������ƣ�my_xtodz
% �������ܣ��غ�����ģ�� ��׼λ�����������λ��
% ���룺xstd      :��ȡλ��    
%      L      :�������λ��
%      u      :������ν�
% �����varargout{1}:z��λ��
%     ��varargout{2}ʵ��λ�����
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
    [~,qD]=my_xtodz(0,L,u);%�������λ�á���������ȡ�غ�ϵ��
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