function dR=my_getdR(R,L,u)
% 函数名称：my_getdR
% 函数功能：相对导航滤波状态转移矩阵
% 输入：R      :当前子惯在主惯坐标系位置
%       L      :标称位置     
%       u      :机翼变形角
% 输出：dR:位置误差量测值
Rf=my_getdRf(L,u,2);%2 表示精确计算
dR=R-L-Rf;