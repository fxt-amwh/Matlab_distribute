function anti=my_getIIlist(dim,T)
% 函数名称：my_getIIlist
% 函数功能：通过原函数在采样周期内的值获得被积函数值
% 输入：dim:当前原函数在单次采样时间内序列
%      T  : 采样时间
% 输出：anti:采样周期内的原函数平均值序列
n=length(dim);
anti=zeros(1,n);
for i=2:n
    anti(i-1)=(dim(i)-dim(i-1))/T;
end