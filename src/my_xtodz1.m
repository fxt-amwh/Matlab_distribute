function [x,dz]=my_xtodz1(xstd,L,u)
% 函数名称：my_xtodz1
% 函数功能：挠曲模型
% 输入：L      :标称位置     
%       u      :机翼变形角
% 输出：Rf:子惯位移
uk=xstd/L*u;
x=xstd.*cos(uk).^2;
dz=xstd.*sin(uk).*cos(uk);