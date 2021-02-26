function INSS=my_invRIpackN(MINS,Smove)
% 函数名称：my_invRIpackN
% 函数功能：反解算相对导航子惯无噪声数据解封 多子不确定数量
% 输入：MINS      :主惯
%       Smove     :子惯运动元组 元组数量对应子惯数量
% 输出：INSS      :子惯结构体元组
[M,N]=size(Smove);%预留升级面阵SINS接口 M=1为线SINS组
INSS=cell(M,N);
for m=1:M
    for n=1:N
        INSS{m,n}=my_invRIpack(MINS,Smove{m,n});
    end
end
