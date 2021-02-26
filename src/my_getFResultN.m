function [INSSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,len,selcCell,select)
% 函数名称：my_SINSgetFResultN
% 函数功能：由子惯输出进行相对导航滤波
% 输入：MINS      :主惯
%       SINS      :子惯
%       KFinit    :子惯滤波器初始值
%       atterr    :子惯安装误差角估计初始值
%       len       :滤波长度       
%       selcCell  :select=1，selcCell表示 SINSR 子惯纯相对导航结果
%                  select=2，selcCell表示 SERR  子惯误差初始值
%       select    :选择子惯误差来源
% 输出：INSSFR:相对导航结果
fprintf('滤波...%5.0f %%',0);
[M,N]=size(SINS);%预留升级面阵SINS接口 M=1为线SINS组
INSSFR=cell(M,N);
Filter=cell(M,N);
if select==1%直接使用相对导航的误差
    SINSR=selcCell;
    for m=1:M
        for n=1:N
        [INSSFR{m,n},Filter{m,n}]=my_getFResult...
            (MINS,SINS{m,n},KFinit{m,n},atterr{m,n},len,...
        SINSR{m,n}.ws_m_addnoise,SINSR{m,n}.fs_m_addnoise);
        end
    end
else%重新注入误差
    SERR=selcCell;
    for m=1:M
        for n=1:N
        [INSSFR{m,n},Filter{m,n}]=my_getFResult...
            (MINS,SINS{m,n},KFinit{m,n},atterr{m,n},len,SERR{m,n});
        end
    end
end
fprintf('\n滤波完成！\n'); 

