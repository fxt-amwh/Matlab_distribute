function y_next=my_LK(my_LKfunc,y_now,u,T,k)
% 函数名称：my_LK
% 函数功能：龙格库塔法求解微分方程数值解
% 输入：my_LKfunc:微分方程函数句柄
%      y_now  : 当前微分项值
%      u      : 函数除微分项其余参数,以元组给入
%       k=4时，u{1}存放当前输入时刻对应的其余参数值
%              u{2}存放中间时刻对应的其余参数值
%              u{3}存放输出时刻对应的其余参数值
%      T      : 采样时间
%      k      : 阶数
% 输出：anti:采样周期内的原函数平均值序列
if k==4
    k1=T*my_LKfunc(y_now,u{1});
    k2=T*my_LKfunc(y_now+0.5*k1,u{2});
    k3=T*my_LKfunc(y_now+0.5*k2,u{2});
    k4=T*my_LKfunc(y_now+k3,u{3});
    y_next=y_now+(k1+2*k2+2*k3+k4)/6;
end
