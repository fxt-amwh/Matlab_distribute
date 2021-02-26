function anti=my_getII(dim,T)
% 函数名称：my_getII
% 函数功能：通过函数在采样周期内的积分获得被积函数值
% 输入：dim:当前被积函数在单次采样时间内的积分序列值
%      T  : 采样时间
% 输出：anti:采样周期内的被积函数平均值序列
if length(dim)==length(T)
    anti=dim./T;
else
    anti=dim/T(1);
end
