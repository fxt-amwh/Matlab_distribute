function k=my_Rtok(R,R0)
% 函数名称：my_getdR
% 函数功能：由R求出挠曲对R的单位比例
% 输入：R      :当前子惯在主惯坐标系位置
%     :R0      :当前子惯标称位置
% 输出：k:位置误差量测值
uf=my_Rtouf(R,R0);
k=[0;uf(2)/sqrt(R0'*R0);0];
% my_getdRf(L,u,flag)


