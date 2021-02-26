%文件名称：test7
%作者信息：冯鑫涛
%功能描述：姿态匹配滤波 有反馈
%版本时间：2021/2/03 19:33
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
SERR.eb=[1;1;1]*dph;
SERR.web=[1;1;1]*dpsh;
SERR.db=[200;200;200]*ug;
SERR.wdb=[200;200;200]*ugpsHz;
%% 速度匹配 无挠曲
close all
R=[2;0;0];
SERR.eb=[1;1;1]*dph;
SERR.web=[1;1;1]*dpsh;
SERR.db=[200;200;200]*ug;
SERR.wdb=[200;200;200]*ugpsHz;
atterr0=[0;0;0]*arcdeg;
SINS_Ret_VLF=zeros(floor(len/2),22);kk=1;t=0;
qnb_vm=a2qua(att0);
vn_vm=vn0;
posMINS_vm=pos0;
qnb_vs=a2qua(att0+atterr0);
vn_vs=vn0;
posSINS_vs=pos0;
cl = cos(posSINS_vs(1,1)); Re = 6378137;
posSINS_vs(2,1)=posSINS_vs(2,1)-R(1)/(cl*Re);
Filter.X=zeros(15,len/2);
Cmserr_1=a2mat(att0);%注意中间处理旋转顺序yxz 不准确的初始角
qns=m2qua(Cmserr_1);%初始姿态 不准确的初始角
fs0=[0;0;0];
KFinit.Qk = diag([ SERR.web;0.001*SERR.web])^2*nts;
KFinit.rk = [0.01;0.01;0.01]*arcdeg;  
KFinit.Rk = diag(KFinit.rk)^2;
KFinit.P0 = diag([[0.01;0.01;0.01]*arcdeg; SERR.eb;[0.1;0.1;0.1]*arcdeg;[0.1;0.1;0.1]*arcdeg;[0.1;0.1;0.1]*arcdeg/nts;])^2;
eth0 = earth(posMINS_vm, vn_vm);  % 地球相关参数计算
K=[0.01;0.01;0.01]*arcdeg;
kfft=my_kfftA15(eth0,q2mat(qnb_vm),q2mat(qns),K,nts);
kf = kfinit(KFinit.Qk, KFinit.Rk, KFinit.P0,kfft.phi,kfft.H);  % kf滤波器初始化
% SERR.eb=[0;0;0]*dph;
% SERR.web=[0;0;0]*dpsh;
% SERR.db=[0;0;0]*ug;
% SERR.wdb=[0;0;0]*ugpsHz;
flag=2;%                                                                    选择有无挠曲
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [ws1,vs1]=imuadderr(SINS{1,flag}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,flag}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL,ethm]=...
        my_Vrelinsupdate...
        (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
    
    kfft=my_kfftA15(ethm,q2mat(qnb_vm),q2mat(qnb_vs),K,nts);
    kf.Phikk_1=kfft.phi;
    kf.Gammak=kfft.Gammak;
    
    ZV_arg=q2mat(qnb_vm)*q2mat(qnb_vs)';% 构造量测
    ZV=[ZV_arg(2,3);ZV_arg(3,1);ZV_arg(1,2)];
    
    kf = kfupdate(kf,ZV,'B');%卡尔曼滤波

%     qnb_vs = qdelphi(qnb_vs,kf.Xk(4:6));  kf.Xk(4:6) = 0;  % ·反馈
%     vn_vs = vn_vs - kf.Xk(1:3);  kf.Xk(1:3) = 0;
    
    Filter.X(:,k)=kf.Xk;
    Filter.X(1:3,k)=q2att(rv2q(Filter.X(1:3,k)));
    Filter.X(7:9,k)=q2att(rv2q(Filter.X(7:9,k)));
    Filter.X(10:12,k)=q2att(rv2q(Filter.X(10:12,k)));
    SINS_Ret_VLF(kk,:)=...
        [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
    kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
Filter.XT=Filter.X';
%% 绘图
% msplot(611, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% msplot(612, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% msplot(613, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,19:21), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% msplot(614, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15)-(SINS_Ret_VLF(:,4:6)), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% msplot(615, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15)-((SINS_Ret_VLF(:,4:6))+SINS_Ret_VLF(:,19:21)), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% msplot(616, SINS_Ret_VLF(:,22), (SINS_Ret_VLF(:,13:15))-((SINS_Ret_VLF(:,4:6))+SINS_Ret_VLF(:,19:21)), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% msplot(211, SINS_Ret_VLF(:,22), (SINS_Ret_VLF(:,1:3)-SINS_Ret_VLF(:,10:12))/arcdeg, 'Att / ( ° )'); legend('\it\theta', '\it\gamma', '\it\psi')
% msplot(212, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)-SINS_Ret_VLF(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%对比滤波估计速度误差
msplot(311, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)+SINS_Ret_VLF(:,19:21), '期望值 Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')%无误差纯惯主惯导结果
msplot(312, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15), '滤波结果 Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(313, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)+SINS_Ret_VLF(:,19:21)-SINS_Ret_VLF(:,13:15), '误差 Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
% 对比滤波估计姿态误差
msplot(611, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,1:3)/arcdeg,'Att / ( ° )'); legend('\it\theta', '\it\gamma', '\it\psi')%无误差纯惯主惯导结果
msplot(612, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,10:12)/arcdeg,'Att / (° )'); legend('\it\theta', '\it\gamma', '\it\psi')%无误差纯惯主惯导结果
msplot(613, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,1:3)/arcdeg-SINS_Ret_VLF(:,10:12)/arcdeg-Filter.XT(:,1:3)/arcdeg+Filter.XT(:,7:9)/arcdeg,'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
ylim([-10,10])
msplot(614, SINS_Ret_VLF(:,22), Filter.XT(:,1:3)/arcdeg,'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(615, SINS_Ret_VLF(:,22), Filter.XT(:,7:9)/arcdeg,'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(616, SINS_Ret_VLF(:,22), Filter.XT(:,10:12)/arcdeg,'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
%% 
msplot(111, SINS_Ret_VLF(:,22), (Filter.XT(:,1:3)-Filter.XT(:,7:9)-Filter.XT(:,10:12))/arcdeg,'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
err_all=(Filter.XT(:,1:3)-Filter.XT(:,7:9)-Filter.XT(:,10:12)+atterr0')/arcdeg;
figure
subplot(311)
plot(SINS_Ret_VLF(:,22),err_all(:,1),'LineWidth',2)
xlabel('时间/s');ylabel('俯仰角/°');title('主子惯导相对姿态误差');grid on;
subplot(312)
plot(SINS_Ret_VLF(:,22),err_all(:,2),'LineWidth',2)
xlabel('时间/s');ylabel('滚转角/°');title('主子惯导相对姿态误差');grid on;
subplot(313)
plot(SINS_Ret_VLF(:,22),err_all(:,3),'LineWidth',2)
xlabel('时间/s');ylabel('航向角/°');title('主子惯导相对姿态误差');grid on;
%% 统计MSE
beg_T=200;
end_T=500;
beg_index=beg_T/nts;
end_index=end_T/nts;
%diag(resultR(:,sam.begin:sam.end)*resultR(:,sam.begin:sam.end)')/length(resultR(:,sam.begin:sam.end))/simLen^2
err=err_all(beg_index:end_index,:)*arcdeg;
MSE = 1000*sqrt(diag(err'*err/length(err)));
fprintf("MSE:(mrad)");
disp(MSE');