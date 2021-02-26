function atterr=my_getAtterr(aList)
% 函数名称：KFinitN
% 函数功能：面阵卡尔曼滤波器
% 输入：SERR      :子惯误差元组
%       SINS      :子惯信息元组  
%       nts = 2*ts;子样数和采样时间
% 输出：atterr:相对导航结果
[~,N,M]=size(aList);%预留升级面阵SINS接口 M=1为线SINS组
atterr=cell(M,N);
for m=1:M
    for n=1:N
        att=aList(:,n,m);
        atterr{m,n}=att;
    end
end
