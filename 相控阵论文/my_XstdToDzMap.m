function [ddx,ddz]=my_XstdToDzMap(xstd,L,u,varargin)
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
[N,M]=size(xstd);
ddx=zeros(N,M);
ddz=zeros(N,M);
for n=1:N
    for m=1:M
        [dz,dx]=my_XstdToDz(xstd(n,m),L,u,err,flag);
        ddx(n,m)=dx;ddz(n,m)=dz;
    end
end