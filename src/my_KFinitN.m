function KFinitCell=my_KFinitN(SERR,SINS,nts)
% 函数名称：KFinitN
% 函数功能：面阵卡尔曼滤波器
% 输入：SERR      :子惯误差元组
%       SINS      :子惯信息元组  
%       nts = 2*ts;子样数和采样时间
% 输出：KFinitCell:滤波器元组
[M,N]=size(SERR);%预留升级面阵SINS接口 M=1为线SINS组
KFinitCell=cell(M,N);
arcdeg = pi/180;
for m=1:M
    for n=1:N
        KFinit.Qk = diag([SERR{m,n}.web; SERR{m,n}.wdb;])^2*nts;
        KFinit.rk = [0.001;0.001;0.001]*sqrt(SINS{m,n}.R0'*SINS{m,n}.R0);  
        KFinit.Rk = diag(KFinit.rk)^2;
        KFinit.P0 = diag([[0.1;0.1;0.1]*arcdeg; [20;20;20]; [0.002;0.002;0.002]*sqrt(SINS{m,n}.R0'*SINS{m,n}.R0);
        SERR{m,n}.eb; SERR{m,n}.db])^2;
        KFinitCell{m,n}=KFinit;
    end
end