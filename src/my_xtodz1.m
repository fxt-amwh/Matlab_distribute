function [x,dz]=my_xtodz1(xstd,L,u)
% �������ƣ�my_xtodz1
% �������ܣ�����ģ��
% ���룺L      :���λ��     
%       u      :������ν�
% �����Rf:�ӹ�λ��
uk=xstd/L*u;
x=xstd.*cos(uk).^2;
dz=xstd.*sin(uk).*cos(uk);