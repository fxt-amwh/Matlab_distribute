%文件名称：test2
%作者信息：冯鑫涛
%功能描述：主子惯导分别进行导航得到速度
%版本时间：2021/1/27 20:23
clc
clear
close all
load MandS_R2_u5.mat
gvar;    % 加载全局变量
g2d=180/pi;
d2g=1/g2d;
att0=MandS.TR.att0;vn0=MandS.TR.vn0;pos0=MandS.TR.pos0;
ts=MandS.MINS.ts;
SINS=MandS.SINS;MINS=MandS.MINS;
TR=MandS.TR;
nn=2;nts=nn*ts;
len=length(SINS{1,2}.wis);
R=[2;0;0];
%% 主惯纯惯导解算
MINS_Ret=zeros(floor(len/2),10);kk=1;t=0;
MERR.eb=[0;0;0]*dph;
MERR.web=[0;0;0]*dpsh;
MERR.db=[0;0;0]*ug;
MERR.wdb=[0;0;0]*ugpsHz;
qnb=a2qua(att0);
vn=vn0;
posMINS=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [qnb,vn,posMINS]=my_insupdate(qnb,vn,posMINS,wm1,vm1,ts);
    MINS_Ret(kk,:)=[q2att(qnb);vn;posMINS;t];kk=kk+1;
    
    if mod(t,100)<nts,disp(fix(t));end
end
%% 子惯刚体纯惯导解算
SINS_Ret_L=zeros(floor(len/2),10);kk=1;t=0;
SERR.eb=[0;0;0]*dph;
SERR.web=[0;0;0]*dpsh;
SERR.db=[0;0;0]*ug;
SERR.wdb=[0;0;0]*ugpsHz;
qnb_L=a2qua(att0);
vn_L=vn0;
posSINS_L=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(SINS{1,1}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,1}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_L,vn_L,posSINS_L]=my_insupdate(qnb_L,vn_L,posSINS_L,wm1,vm1,ts);
    SINS_Ret_L(kk,:)=[q2att(qnb_L);vn_L;posSINS_L;t];kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
%% 子惯挠曲纯惯导解算
SINS_Ret_u=zeros(floor(len/2),10);kk=1;t=0;
qnb_u=a2qua(att0);
vn_u=vn0;
posSINS_u=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(SINS{1,2}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,2}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_u,vn_u,posSINS_u]=my_insupdate(qnb_u,vn_u,posSINS_u,wm1,vm1,ts);
    SINS_Ret_u(kk,:)=[q2att(qnb_u);vn_L;posSINS_u;t];kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
%% 绘图
delpos=deltapos(TR.pos);
delposMINS=deltapos(MINS_Ret(:,7:9));
delposSINS_L=deltapos(SINS_Ret_L(:,7:9));
delposSINS_u=deltapos(SINS_Ret_u(:,7:9));
%主惯轨迹与子惯刚体误差
%% 需要使用的速度
msplot(211, MINS_Ret(:,10), MINS_Ret(:,1:3)/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, MINS_Ret(:,10), MINS_Ret(:,4:6), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(211, SINS_Ret_L(:,10), SINS_Ret_L(:,1:3)/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_L(:,10), SINS_Ret_L(:,4:6), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(211, SINS_Ret_u(:,10), SINS_Ret_u(:,1:3)/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_u(:,10), SINS_Ret_u(:,4:6), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')

msplot(211, SINS_Ret_u(:,10), (SINS_Ret_L(:,1:3)-MINS_Ret(:,1:3))/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_u(:,10), SINS_Ret_L(:,4:6)-MINS_Ret(:,4:6), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%% 速度匹配 无挠曲
SINS_Ret_VL=zeros(floor(len/2),22);kk=1;t=0;
qnb_vm=a2qua(att0);
vn_vm=vn0;
posMINS_vm=pos0;
qnb_vs=a2qua(att0);
vn_vs=vn0;
posSINS_vs=pos0;
cl = cos(posSINS_vs(1,1)); Re = 6378137;
posSINS_vs(2,1)=posSINS_vs(2,1)+R(1)/(cl*Re);
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [ws1,vs1]=imuadderr(SINS{1,1}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,1}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL]=...
        my_Vrelinsupdate...
        (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
    SINS_Ret_VL(kk,:)=...
        [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
    kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
% 绘图
msplot(211, SINS_Ret_VL(:,22), (SINS_Ret_VL(:,1:3)-SINS_Ret_VL(:,1:3))/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_VL(:,22), SINS_Ret_VL(:,4:6)+SINS_Ret_VL(:,19:21)-SINS_Ret_VL(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%% 速度匹配 有挠曲
SINS_Ret_Vu=zeros(floor(len/2),22);kk=1;t=0;
qnb_vm=a2qua(att0);
vn_vm=vn0;
posMINS_vm=pos0;
qnb_vs=a2qua(att0);
vn_vs=vn0;
posSINS_vs=pos0;
cl = cos(posSINS_vs(1,1)); Re = 6378137;
posSINS_vs(2,1)=posSINS_vs(2,1)+R(1)/(cl*Re);
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [ws1,vs1]=imuadderr(SINS{1,2}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,2}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL]=...
        my_Vrelinsupdate...
        (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
    SINS_Ret_Vu(kk,:)=...
        [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
    kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
% 绘图
msplot(211, SINS_Ret_Vu(:,22), (SINS_Ret_Vu(:,1:3)-SINS_Ret_Vu(:,1:3))/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_Vu(:,22), SINS_Ret_Vu(:,4:6)+SINS_Ret_Vu(:,19:21)-SINS_Ret_Vu(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')