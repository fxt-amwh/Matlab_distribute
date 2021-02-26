function INSSR=my_SINSgetResultN(MINS,SINS,len,SERR,atterr)
% 函数名称：my_SINSgetResultN
% 函数功能：由子惯输出进行相对导航
% 输入：MINS   :主惯导
%       SINS   :子惯导数据元组 M*N 
%       len    :解算序列长度的一半 1*1
%       SERR   :子惯导误差元组 M*N 
%       atterr :子惯导位置误差元组 M*N 
% 输出：INSSR:相对导航结果 面接口
%% 相对导航
fprintf('纯子惯相对导航...%5.0f %%',0); 
[M,N]=size(SINS);%预留升级面阵SINS接口 M=1为线SINS组
INSSR=cell(M,N);
for m=1:M
    for n=1:N
        INSSR{m,n}=my_SINSgetResult(MINS,SINS{m,n},len,SERR{m,n},atterr{m,n});
    end
end
fprintf('\n纯惯相对导航完成！\n'); 