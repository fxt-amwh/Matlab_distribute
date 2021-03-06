%% 一主 + 一子 重要
%文件名称：test6
%作者信息：冯鑫涛
%功能描述：相对导航滤波 速度匹配对 有反馈
%版本时间：2021/3/15 19:15
% 1、根据设计轨迹运动反解算主惯传感器输出
% 2、根据设计机翼运动反解算子惯传感器数据
% 3、子惯数据注入误差
% 4、纯惯性相对导航
% 5、滤波
% 6、误差计算
%% 主惯运动数据设计仿真
clc
clear
close all
gvar;    % 加载全局变量
ts = 0.01;%采样时间
att0 = [0;0;0]*arcdeg; vn0 = [0;100;0]; pos0 = [[34;108]*arcdeg;100];
%     俯仰角速率 横滚角速率 方位角速率 纵向加速度 持续时间 （环形轨迹设计）
wat = [  0,         0,          0,         0,         10       %静止
        60,         0,          0,         0,         60       %加速
        -60,        0,          0,         0,         60       %匀速
        0,          0,          0,         0,         50       %匀速
        0,          0,          60,       0,         360      %匀速转弯
        0,          0,          0,         0,         50       %匀速
        -60,        0,          0,         0,         60       %匀速
        60,         0,          0,         0,        60       %减速
        0,          0,          0,         0,         10       %减速
        ];    %静止
wat(:,1:3) = wat(:,1:3)*arcdeg/60;  % deg/min->deg/s
fprintf('轨迹仿真中...'); 
[attm, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
fprintf('轨迹仿真完成！\n主惯仿真中...');
[wm, vm] = my_av2imu(attm, vn, pos, ts);
tt = (0:length(attm)-1)'*ts;
%% 主惯信息
wim0=my_getII(wm(1,:)',ts);
wim=[wim0,my_getII(wm',ts)];%wim注意是在主坐标系下的量测
fm=[my_getII(vm(1,:)',ts),my_getII(vm',ts)];
MINS.wim=wim;
MINS.wim0=wim0;
MINS.fm=fm;
fprintf('主惯仿真完成！\n');
%% 子惯运动数据设计
fprintf('子惯仿真...'); 
len=size(wm,1)+1;% len=7200;
um=2*arcdeg;
f=0.01;
u=[zeros(1,len);um*sin(2*pi*f*(0:ts:((len-1)*ts)));zeros(1,len)];
R0=[2;0;0];%理想时 主子间0时刻前初始相对位置
Rf0=[0;0;0];%理想时 子0时刻前偏移
Rf=[-R0(1)*sin(u(2,:)).*sin(u(2,:));zeros(1,len);R0(1)*sin(u(2,:)).*cos(u(2,:))];
R=R0+Rf;
att=2*u;
att_s0=[0;0;0]*arcdeg;%理想时 主子间0时刻前初始相对角度
%% 子惯反解算数据
% atterr0=[1;2;3]*arcdeg;%注意这里旋转顺序zxy
atterr0=[0;0;0]*arcdeg;%注意这里旋转顺序zxy
Cmserr=my_a2mat(atterr0);%注意中间处理旋转顺序yxz
diffRf0=[0;0;0];%Rf的初始变化率
U0=diffRf0+cross(wim0,R0);
SINS1=my_invRI(R,att,atterr0,wim,fm,ts,R0,att_s0,U0);
fprintf('子惯仿真完成！\n纯子惯相对导航...%5.0f %%',0);  
%% 子惯误差设计
SERR1.eb = [1;1;1]*dph; SERR1.web = [1;1;1]*dpsh;   %陀螺常值零偏，角度随机游走
SERR1.db = [200;200;200]*ug; SERR1.wdb = [200;200;200]*ugpsHz;  %加速度计常值偏值，速度随机游走
% SERR1.eb = 0*[1;1;1]*dph; SERR1.web = 0*[1;1;1]*dpsh;   %陀螺常值零偏，角度随机游走
% SERR1.db = 0*[200;200;200]*ug; SERR1.wdb = 0*[200;200;200]*ugpsHz;  %加速度计常值偏值，速度随机游走
%% 相对导航
atterr0_1=[10;5;2]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
SINSR1=my_SINSgetResult(MINS,SINS1,len/2,SERR1,atterr0_1);
fprintf('纯惯相对导航完成！\n滤波...%5.0f %%',0);
%% 滤波
nts = 2*ts; %子样数和采样时间
KFinit1.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit1.rk = [0.01;0.01;0.01];  
KFinit1.Rk = diag(KFinit1.rk)^2;
KFinit1.P0 = diag([[0.1;0.1;0.1]*arcdeg; [0.001;0.001;0.001]; [0.001;0.001;0.001];
         [1;1;1]*dph; [200;200;200]*ug])^2;
[SINSFR1,~]=my_getFResult(MINS,SINS1,KFinit1,atterr0_1,len/2,...
    SINSR1.ws_m_addnoise,SINSR1.fs_m_addnoise);
fprintf('滤波完成！\n'); 
%% 角度真值 zxy转序
atttrue=zeros(3,len);
for i=1:len%真实角度
   [rz,rx,ry]=dcm2angle(SINS1.Cms{i},'zxy');
   atttrue(:,i)=[rx,ry,-rz]; 
end
%% 速度匹配 无挠曲
MERR.eb=[0;0;0]*dph;
MERR.web=[0;0;0]*dpsh;
MERR.db=[0;0;0]*ug;
MERR.wdb=[0;0;0]*ugpsHz;
Flag.EnBack=true;
SINS_Ret_VLF=zeros(floor(len/2),22);kk=1;t=0;
qnb_vm=a2qua(att0);
vn_vm=vn0;
posMINS_vm=pos0;
qnb_vs=a2qua(att0+atterr0_1);
vn_vs=vn0;
posSINS_vs=pos0;
cl = cos(posSINS_vs(1,1)); Re = 6378137;
posSINS_vs(2,1)=posSINS_vs(2,1)-R0(1)/(cl*Re);
Filter.X=zeros(12,len/2);
Cmserr_1=a2mat(att0+atterr0_1);%注意中间处理旋转顺序yxz 不准确的初始角
qns=m2qua(Cmserr_1);%初始姿态 不准确的初始角
fs0=[0;0;0];
KFinit.Qk = diag([SERR1.wdb; SERR1.web;])^2*nts;
KFinit.rk = [0.001;0.001;0.001];  
KFinit.Rk = diag(KFinit.rk)^2;
KFinit.P0 = diag([[0.001;0.001;0.001];[0.1;0.1;0.1]*arcdeg;SERR1.db; SERR1.eb])^2;
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
    [ws1,vs1]=imuadderr(SINS1.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS1.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR1.eb, SERR1.web, SERR1.db, SERR1.wdb, ts);
    [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL,ethm]=...
        my_Vrelinsupdate...
        (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R0,ts);
    
    fs=sum(SINS1.fs(:,(2*k-1):(2*k)),2)/2;
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
figure%理论相对位置与纯解算相对位置对比
subplot(311)
T=ts:ts:(len*ts);
plot(T(2:4:end),R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,1),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('x/m');
title('主子惯导相对位置R（滤波）');
grid on;
subplot(312)
plot(T(2:4:end),R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,2),'--','LineWidth',2)
grid on;
legend('真实值','计算值');
xlabel('时间/s');ylabel('y/m');
title('主子惯导相对位置R（滤波）');
subplot(313)
plot(T(2:4:end),R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,3),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('z/m');
title('主子惯导相对位置R（滤波）');
grid on;
%%
figure%理论角度与纯相对解算角度对比
subplot(311)
plot(T(2:4:end-1),(atttrue(1,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(1,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('俯仰角/°');
title('主子惯导相对姿态（滤波）');

grid on;
subplot(312)
plot(T(2:4:end-1),(atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(2,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('滚转角/°');
title('主子惯导相对姿态（滤波）');

grid on;
subplot(313)
plot(T(2:4:end-1),(atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(3,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('航向角/°');
title('主子惯导相对姿态（滤波）');
grid on;
%% 
figure
subplot(311)
T=ts:ts:(len*ts);
plot(T(2:4:end),SINSFR1.Rall(1:2:end,1)'-R(1,2:4:end),'LineWidth',2)
legend('X');
xlabel('时间/s');ylabel('x/m');
title('主子惯导相对位置误差');
grid on;
subplot(312)
plot(T(2:4:end),SINSFR1.Rall(1:2:end,2)'-R(2,2:4:end),'LineWidth',2)
grid on;
legend('Y');
xlabel('时间/s');ylabel('y/m');
title('主子惯导相对位置误差');
subplot(313)
plot(T(2:4:end),SINSFR1.Rall(1:2:end,3)'-R(3,2:4:end),'LineWidth',2)
legend('Z');
xlabel('时间/s');ylabel('z/m');
title('主子惯导相对位置误差');
grid on;
%%
figure
subplot(311)
plot(T(2:4:end-1),(SINSFR1.attall(1,1:2:end-1)-atttrue(1,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),(Filter.XT(1:2:end,4)+atterr0_1(1,1))'/arcdeg,'LineWidth',2)
legend("相对导航","速度匹配");
xlabel('时间/s');ylabel('俯仰角/°');title('主子惯导相对姿态误差');grid on;
subplot(312)
plot(T(2:4:end-1),(SINSFR1.attall(2,1:2:end-1)-atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),(Filter.XT(1:2:end,5))'/arcdeg,'LineWidth',2)
legend("相对导航","速度匹配");
xlabel('时间/s');ylabel('滚转角/°');title('主子惯导相对姿态误差');grid on;
ylim([-1,1]);
subplot(313)
plot(T(2:4:end-1),(SINSFR1.attall(3,1:2:end-1)-atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),(Filter.XT(1:2:end,5))'/arcdeg,'LineWidth',2)
legend("相对导航","速度匹配");
xlabel('时间/s');ylabel('航向角/°');title('主子惯导相对姿态误差');grid on;
ylim([-1,1]);
%%
figure
subplot(311)
plot(T(2:4:end-1),(SINSFR1.attall(1,1:2:end-1)-atttrue(1,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),-SINS_Ret_VLF(1:2:end,1)/arcdeg+SINS_Ret_VLF(1:2:end,10)/arcdeg,'LineWidth',2)
legend("相对导航","速度匹配");
xlabel('时间/s');ylabel('俯仰角/°');title('主子惯导相对姿态误差');grid on;
subplot(312)
plot(T(2:4:end-1),(SINSFR1.attall(2,1:2:end-1)-atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINS_Ret_VLF(1:2:end,2)/arcdeg-SINS_Ret_VLF(1:2:end,11)/arcdeg,'LineWidth',2)
legend("相对导航","速度匹配");
xlabel('时间/s');ylabel('滚转角/°');title('主子惯导相对姿态误差');grid on;
ylim([-1,1]);
subplot(313)
plot(T(2:4:end-1),(SINSFR1.attall(3,1:2:end-1)-atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINS_Ret_VLF(1:2:end,3)/arcdeg-SINS_Ret_VLF(1:2:end,12)/arcdeg,'LineWidth',2)
legend("相对导航","速度匹配");
xlabel('时间/s');ylabel('航向角/°');title('主子惯导相对姿态误差');grid on;
ylim([-1,1]);
%% 统计MSE
beg_T=200;
end_T=500;
beg_index=beg_T/nts;
end_index=end_T/nts;
%diag(resultR(:,sam.begin:sam.end)*resultR(:,sam.begin:sam.end)')/length(resultR(:,sam.begin:sam.end))/simLen^2
errall1=SINSFR1.attall(:,1:end)-atttrue(:,1:2:end);
err1=errall1(:,beg_index:end_index);
MSE1 = 1000*sqrt(diag(err1*err1'/length(err1)));
errall2=SINS_Ret_VLF(1:end,1:3)-SINS_Ret_VLF(1:end,10:12)+atttrue(:,1:2:end)';
errall2=errall2';
if  Flag.EnBack
    err2=errall2(:,beg_index:end_index);
    err2(3,err2(3,:)<-pi)=0.1*pi/180;err2(3,err2(3,:)>pi)=0.1*pi/180;
    MSE2 = 1000*sqrt(diag(err2*err2'/length(err2)));
else    
    err2=Filter.XT(beg_index:end_index,4:6)+atterr0_1';
    MSE2 = 1000*sqrt(diag(err2'*err2/length(err2)));
end
MSE1=[MSE1;sqrt(MSE1'*MSE1)];
MSE2=[MSE2;sqrt(MSE2'*MSE2)];
fprintf("MSE相对导航:(mrad)");
disp(MSE1');
fprintf("MSE速度匹配:(mrad)");
disp(MSE2');