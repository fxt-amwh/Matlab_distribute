function INSS=my_invRIpack(MINS,Smove)
% 函数名称：my_invRIpack
% 函数功能：反解算相对导航子惯无噪声数据解封
% 输入：MINS      :主惯
%       Smove     :子惯运动
% 输出：INSS      :子惯结构体
INSS=my_invRI(Smove.R,Smove.att,Smove.atterr0,...
    MINS.wim,MINS.fm,MINS.ts,Smove.R0,Smove.att_s0,Smove.U0);
