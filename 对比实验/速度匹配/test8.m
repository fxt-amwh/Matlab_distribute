%文件名称：test4 重要 挠曲对比
%作者信息：冯鑫涛
%功能描述：速度匹配滤波 有反馈
%版本时间：2021/1/28 19:43
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
MERR.eb=[0;0;0]*dph;
MERR.web=[0;0;0]*dpsh;
MERR.db=[0;0;0]*ug;
MERR.wdb=[0;0;0]*ugpsHz;
%% 子惯刚体纯惯导解算
% SERR.eb=[0;0;0]*dph;
% SERR.web=[0;0;0]*dpsh;
% SERR.db=[0;0;0]*ug;
% SERR.wdb=[0;0;0]*ugpsHz;
SERR.eb=[1;1;1]*dph;
SERR.web=[1;1;1]*dpsh;
SERR.db=[200;200;200]*ug;
SERR.wdb=[200;200;200]*ugpsHz;
% %% 速度匹配 无挠曲
% SINS_Ret_VL=zeros(floor(len/2),22);kk=1;t=0;
% qnb_vm=a2qua(att0);
% vn_vm=vn0;
% posMINS_vm=pos0;
% qnb_vs=a2qua(att0);
% vn_vs=vn0;
% posSINS_vs=pos0;
% cl = cos(posSINS_vs(1,1)); Re = 6378137;
% posSINS_vs(2,1)=posSINS_vs(2,1)+R(1)/(cl*Re);
% for k=1:len/2
%     t=t+nts;
%     [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
%         MINS.fm(:,(2*k-1):(2*k))'*ts, ...
%         MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
%     [ws1,vs1]=imuadderr(SINS{1,1}.wis(:,(2*k-1):(2*k))'*ts, ...
%         SINS{1,1}.fs(:,(2*k-1):(2*k))'*ts, ...
%         SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
%     [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL]=...
%         my_Vrelinsupdate...
%         (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
%     SINS_Ret_VL(kk,:)=...
%         [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
%     kk=kk+1;
%     if mod(t,100)<nts,disp(fix(t));end
% end
% %% 绘图
% msplot(211, SINS_Ret_VL(:,22), (SINS_Ret_VL(:,1:3)-SINS_Ret_VL(:,10:12))/arcdeg, 'Att / ( ° )'); legend('\it\theta', '\it\gamma', '\it\psi')
% msplot(212, SINS_Ret_VL(:,22), SINS_Ret_VL(:,4:6)+SINS_Ret_VL(:,19:21)-SINS_Ret_VL(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% msplot(211, SINS_Ret_VL(:,22), (SINS_Ret_VL(:,1:3)-SINS_Ret_VL(:,10:12))/arcdeg, 'Att / ( ° )'); legend('\it\theta', '\it\gamma', '\it\psi')
% msplot(212, SINS_Ret_VL(:,22), SINS_Ret_VL(:,4:6)-SINS_Ret_VL(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% %% 速度匹配 有挠曲
% SINS_Ret_Vu=zeros(floor(len/2),22);kk=1;t=0;
% qnb_vm=a2qua(att0);
% vn_vm=vn0;
% posMINS_vm=pos0;
% qnb_vs=a2qua(att0);
% vn_vs=vn0;
% posSINS_vs=pos0;
% cl = cos(posSINS_vs(1,1)); Re = 6378137;
% posSINS_vs(2,1)=posSINS_vs(2,1)+R(1)/(cl*Re);
% for k=1:len/2
%     t=t+nts;
%     [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
%         MINS.fm(:,(2*k-1):(2*k))'*ts, ...
%         MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
%     [ws1,vs1]=imuadderr(SINS{1,2}.wis(:,(2*k-1):(2*k))'*ts, ...
%         SINS{1,2}.fs(:,(2*k-1):(2*k))'*ts, ...
%         SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
%     [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL]=...
%         my_Vrelinsupdate...
%         (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
%     SINS_Ret_Vu(kk,:)=...
%         [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
%     kk=kk+1;
%     if mod(t,100)<nts,disp(fix(t));end
% end
% %% 绘图
% msplot(211, SINS_Ret_Vu(:,22), (SINS_Ret_Vu(:,1:3)-SINS_Ret_Vu(:,10:12))/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
% msplot(212, SINS_Ret_Vu(:,22), SINS_Ret_Vu(:,4:6)+SINS_Ret_Vu(:,19:21)-SINS_Ret_Vu(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%% 速度匹配 无挠曲
Flag.EnBack=false;
close all
R=[2;0;0];
SERR.eb=[1;1;1]*dph;
SERR.web=[1;1;1]*dpsh;
SERR.db=[200;200;200]*ug;
SERR.wdb=[200;200;200]*ugpsHz;
atterr0=[0;0;10]*arcdeg;
SINS_Ret_VLF=zeros(floor(len/2),22);kk=1;t=0;
qnb_vm=a2qua(att0);
vn_vm=vn0;
posMINS_vm=pos0;
qnb_vs=a2qua(att0+atterr0);
vn_vs=vn0;
posSINS_vs=pos0;
cl = cos(posSINS_vs(1,1)); Re = 6378137;
posSINS_vs(2,1)=posSINS_vs(2,1)-R(1)/(cl*Re);
Filter.X=zeros(12,len/2);
Cmserr_1=a2mat(att0+atterr0);%注意中间处理旋转顺序yxz 不准确的初始角
qns=m2qua(Cmserr_1);%初始姿态 不准确的初始角
fs0=[0;0;0];
KFinit.Qk = diag([SERR.wdb; SERR.web;])^2*nts;
KFinit.rk = [0.001;0.001;0.001];  
KFinit.Rk = diag(KFinit.rk)^2;
KFinit.P0 = diag([[0.001;0.001;0.001];[0.1;0.1;0.1]*arcdeg;SERR.db; SERR.eb])^2;
eth0 = earth(posMINS_vm, vn_vm);  % 地球相关参数计算
kfft=my_kfftV12(eth0,q2mat(qns),fs0,nts);
kf = kfinit(KFinit.Qk, KFinit.Rk, KFinit.P0,kfft.phi,kfft.H);  % kf滤波器初始化
% SERR.eb=[0;0;0]*dph;
% SERR.web=[0;0;0]*dpsh;
% SERR.db=[0;0;0]*ug;
% SERR.wdb=[0;0;0]*ugpsHz;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [ws1,vs1]=imuadderr(SINS{1,2}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,2}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL,ethm]=...
        my_Vrelinsupdate...
        (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
    
    fs=sum(SINS{1,2}.fs(:,(2*k-1):(2*k)),2)/2;
    kfft=my_kfftV12(ethm,q2mat(qnb_vs),fs,nts);
    kf.Phikk_1=kfft.phi;
    kf.Gammak=kfft.Gammak;
    
    ZV=vn_vs-(vn_vm+VL);% 构造量测
    
    kf = kfupdate(kf,ZV,'B');%卡尔曼滤波

    if Flag.EnBack
        qnb_vs = qdelphi(qnb_vs,kf.Xk(4:6));  kf.Xk(4:6) = 0;  % ·反馈
        vn_vs = vn_vs - kf.Xk(1:3);  kf.Xk(1:3) = 0;
    end
    
    Filter.X(:,k)=kf.Xk;
    Filter.X(4:6,k)=q2att(rv2q(Filter.X(4:6,k)));
    SINS_Ret_VLF(kk,:)=...
        [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
    kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
Filter.XT=Filter.X';
%% 绘图
msplot(611, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U');title('主惯导速度');
msplot(612, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U');title('子惯导速度');
msplot(613, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,19:21), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U');title('杆臂速度');
msplot(614, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15)-(SINS_Ret_VLF(:,4:6)), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(615, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15)-((SINS_Ret_VLF(:,4:6))+SINS_Ret_VLF(:,19:21)), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U');title('导航解算主子惯导速度差');
msplot(616, SINS_Ret_VLF(:,22), (SINS_Ret_VLF(:,13:15)-Filter.XT(:,1:3))-((SINS_Ret_VLF(:,4:6))+SINS_Ret_VLF(:,19:21)), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U');title('导航解算主子惯导速度差');
% msplot(211, SINS_Ret_VLF(:,22), (SINS_Ret_VLF(:,1:3)-SINS_Ret_VLF(:,10:12))/arcdeg, 'Att / ( ° )'); legend('\it\theta', '\it\gamma', '\it\psi')
% msplot(212, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)-SINS_Ret_VLF(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%对比滤波估计速度误差
msplot(311, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)+SINS_Ret_VLF(:,19:21), '子惯速度期望值 Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
if Flag.EnBack
    msplot(312, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15), '子惯速度纯惯结果 Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
else
    msplot(312, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15)-Filter.XT(:,1:3), '子惯速度滤波结果 Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
end
msplot(313, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)+SINS_Ret_VLF(:,19:21)-SINS_Ret_VLF(:,13:15)+Filter.XT(:,1:3), '误差 Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% 对比滤波估计姿态误差
msplot(411, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,1:3)/arcdeg,'期望姿态Att / ( ° )'); legend('\it\theta', '\it\gamma', '\it\psi')%无误差纯惯主惯导结果
if Flag.EnBack
    msplot(412, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,10:12)/arcdeg,'子惯姿态滤波结果Att / (° )'); legend('\it\theta', '\it\gamma', '\it\psi')%无误差纯惯主惯导结果
else
    msplot(412, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,10:12)/arcdeg,'子惯纯惯姿态Att / (° )'); legend('\it\theta', '\it\gamma', '\it\psi')%无误差纯惯主惯导结果
end
msplot(413, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,1:3)/arcdeg-SINS_Ret_VLF(:,10:12)/arcdeg,'1-2Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
ylim([-0.5,0.5])
msplot(414, SINS_Ret_VLF(:,22), Filter.XT(:,4:6)/arcdeg,'滤波姿态误差估计Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
%%
if Flag.EnBack
    msplot(111, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,1:3)/arcdeg-SINS_Ret_VLF(:,10:12)/arcdeg+SINS{1,2}.atttrue(:,1:2:end)'/arcdeg,'滤波姿态误差估计Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
    msplot(111, SINS_Ret_VLF(:,22), SINS{1,2}.atttrue(:,1:2:end)/arcdeg,'滤波姿态误差估计Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
else
    msplot(111, SINS_Ret_VLF(:,22), Filter.XT(:,4:6)/arcdeg,'滤波姿态误差估计Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(111, SINS_Ret_VLF(:,22), SINS{1,2}.atttrue(:,1:2:end)/arcdeg,'滤波姿态误差估计Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
end
%% 统计MSE
beg_T=200;
end_T=500;
beg_index=beg_T/nts;
end_index=end_T/nts;
%diag(resultR(:,sam.begin:sam.end)*resultR(:,sam.begin:sam.end)')/length(resultR(:,sam.begin:sam.end))/simLen^2
if Flag.EnBack
    errL=SINS_Ret_VLF(:,1:3)-SINS_Ret_VLF(:,10:12)+SINS{1,2}.atttrue(:,1:2:end)';
    err=errL(beg_index:end_index,:);
    err(err(:,3)<-pi/2,3)=0;
else
    err=Filter.XT(beg_index:end_index,4:6)+atterr0';
end
MSE = 1000*sqrt(diag(err'*err/length(err)));
fprintf("MSE:(mrad)");
disp(MSE');