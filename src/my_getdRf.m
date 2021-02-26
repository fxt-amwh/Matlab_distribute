function Rf=my_getdRf(L,u,flag)
% 函数名称：my_getdRf
% 函数功能：相对导航滤波状态转移矩阵
% 输入：L      :标称位置     
%       u      :机翼变形角
% 输出：Rf:子惯位移
if nargin<3
    flag=1;
end
ud2=u/2;
switch flag
    case 1%小角度近似
        Rf=[-L(1)*(ud2(2)^2)+L(2)*ud2(3)
            -L(2)*(ud2(3)^2)+L(3)*ud2(1)
            -L(3)*(ud2(1)^2)+L(1)*ud2(2)];
    case 2%精确计算
        Rf=[-L(1)*sin(ud2(2))*sin(ud2(2))+L(2)*sin(ud2(3))*cos(ud2(3))
            -L(2)*sin(ud2(3))*sin(ud2(3))+L(3)*sin(ud2(1))*cos(ud2(1))
            -L(3)*sin(ud2(1))*sin(ud2(1))+L(1)*sin(ud2(2))*cos(ud2(2))];
end