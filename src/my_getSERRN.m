function SERR=my_getSERRN(SERRList)
% 函数名称：KFinitN
% 函数功能：面阵卡尔曼滤波器
% 输入：SERR      :子惯误差元组
%       SINS      :子惯信息元组  
%       nts = 2*ts;子样数和采样时间
% 输出：atterr:相对导航结果
[~,N,M]=size(SERRList.eb);%预留升级面阵SINS接口 M=1为线SINS组
SERR=cell(M,N);
for m=1:M
    for n=1:N
        s.eb = SERRList.eb(:,n,m); s.web = SERRList.web(:,n,m);   %陀螺常值零偏，角度随机游走
        s.db = SERRList.db(:,n,m); s.wdb = SERRList.wdb(:,n,m);  %加速度计常值偏值，速度随机游走
        SERR{m,n}=s;
    end
end
