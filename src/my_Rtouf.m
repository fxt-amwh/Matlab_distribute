function uf=my_Rtouf(R,R0)
% 函数名称：my_getdR
% 函数功能：由R求出挠曲角
% 输入：R      :当前子惯在主惯坐标系位置
%     :R0      :当前子惯标称位置
% 输出：uf:位置误差量测值
Rf=R-R0;
if Rf(3)>0
    uf=[0;asin(sqrt((Rf'*Rf)/(R0'*R0)));0];
else
    uf=[0;-asin(sqrt((Rf'*Rf)/(R0'*R0)));0];
end
