function varargout=my_XstdToDzAndAtt(xstd,L,u,varargin)
% �������ƣ�my_XstdToDzAndAtt
% �������ܣ��غ�����ģ�� ��׼λ�����������λ�� ������
% ���룺xstd      :��ȡλ��    
%      L      :�������λ��
%      u      :������ν�
% �����varargout{1}:z��λ��
%     ��varargout{2}ʵ��λ�����
%     ��varargout{3}i��������
%     ��varargout{4}xstd��ʱ��������
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
    [~,qD]=my_xtodz(0,L,u);%�������λ�á���������ȡ�غ�ϵ�� %ע�������qD
end
fw=@(x)sqrt(1+qD*(4*x.^3-12*L*x.^2+12*L^2*x));%���߻���
xloop=xstd;
dx=1;
i=0;
while abs(dx) > err
    dx=integral(fw,0,xloop)-xstd;
    xloop=xloop-dx;
    i=i+1;
end

difffw=@(x)(qD*(4*x.^3-12*L*x.^2+12*L^2*x));%�Ƕ�΢��


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
if 4==nargout
    varargout{1}=my_xtodz(xloop,L,u,qD);
    varargout{2}=xloop-xstd;
    varargout{3}=i;
    varargout{4}=atan(difffw(xloop));
end
if 5==nargout
    varargout{1}=my_xtodz(xloop,L,u,qD);
    varargout{2}=xloop-xstd;
    varargout{3}=i;
    varargout{4}=atan(difffw(xloop));
    varargout{5}=atan(difffw(xstd));
end
