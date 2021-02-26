%% 一主+2子（中间备份）
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
ts = 0.1;%采样时间
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;105];
%     俯仰角速率 横滚角速率 方位角速率 纵向加速度 持续时间 （环形轨迹设计）
% wat = [ 0,         -60,        120,       0,         30       %静止30
%         0,         60,         -120,      0,         60       %加速90
%         0,         0,          0,         1,         20       %匀速110
%         60,        0,          0,         0,         30        %匀速140
%         0,         0,          0,         0,         10       %匀速150
%         -60,       0,          0,         0,         30        %匀速180
%         0,         -60,        120,       0,         60       %静止240
%         0,         60,         -120,      0,         60       %加速300
%         0,         0,          0,         0,         30       %匀速330
%         0,         0,          60,        -1,        20       %匀速350
%         60,        0,          0,         0,         30        %匀速380
%         0,         0,          0,         0,         10       %匀速390
%         -60,       0,          0,         0,         30        %匀速420
%         0,         0,          0,         0,         60       %匀速480
%         -60,       0,          0,         0,         30        %匀速510
%         0,         0,          0,         0,         10       %匀速520
%         60,        0,          0,         0,         30        %匀速550
%         0,         0,          0,         0,         60       %匀速610
%         0,         0,          60,        0,         20        %匀速630
%         0,         0,          -60,       0,         20        %匀速650
%         -60,       0,          0,         0,         30        %匀速680
%         0,         0,          0,         0,         10       %匀速690
%         60,        0,          0,         0,         30        %匀速720
%         
%         ];    %静止
wat = [  0,         0,          0,         0,         10       %静止10
        60,         0,          0,         0,         60       %加速70
        -60,        0,          0,         0,         60       %匀速130
        0,          0,          0,         0,         50       %匀速180
        0,          0,          60,       0,         360      %匀速转弯540
        0,          0,          0,         0,         50       %匀速590
        -60,        0,          0,         0,         60       %匀速650
        60,         0,          0,         0,        60       %减速710
        0,          0,          0,         0,         10       %减速720
        ];    %静止
wat(:,1:3) = wat(:,1:3)*arcdeg/60;  % ->deg/s
fprintf('轨迹仿真中...'); 
[attm, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
fprintf('轨迹仿真完成！\n主惯仿真中...');
[wm, vm] = my_av2imu(attm, vn, pos, ts);
tt = (0:length(attm)-1)'*ts;
% 轨迹作图
msplot(221, tt, attm/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(222, tt, vn, 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(223, tt, deltapos(pos), '\DeltaPos / m');
      legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')
msplot(224, pos(:,2)/arcdeg, pos(:,1)/arcdeg, '\itL\rm / ( \circ )', '\it\lambda\rm / ( \circ)');
      hold on, plot(pos(1,2)/arcdeg, pos(1,1)/arcdeg, 'ro');
%  惯性传感器信息作图
msplot(121, tt(2:end), wm/ts/arcdeg, '\it\omega^b_{ib}\rm / ( \circ.s^{-1} )');
      legend('\it\omega^b_{ibx}', '\it\omega^b_{iby}', '\it\omega^b_{ibz}');
msplot(122, tt(2:end), vm/ts, '\itf^b\rm_{sf} / ( m.s^{-2} )');
      legend('\itf^b\rm_{sf\itx}', '\itf^b\rm_{sf\ity}', '\itf^b\rm_{sf\itz}');
% 三维轨迹
%%
figure
delpos=deltapos(pos);
plot3(pos(:,1)/arcdeg,pos(:,2)/arcdeg,delpos(:,3));
% axis([0,-8000,]);
grid on;
xlabel('纬度/°'),ylabel('经度/°'),zlabel('高/m');
%% 主惯信息
wim0=my_getII(wm(1,:)',ts);
wim=[wim0,my_getII(wm',ts)];%wim注意是在主坐标系下的量测
fm=[my_getII(vm(1,:)',ts),my_getII(vm',ts)];
MINS.wim=wim;
MINS.wim0=wim0;
MINS.fm=fm;
MINS.ts=ts;
fprintf('主惯仿真完成！\n');
%% 子惯运动数据设计
fprintf('子惯仿真...'); 
Sinf.ts=MINS.ts;
len=size(wm,1)+1;
Sinf.len=size(wm,1)+1;% len=7200;
%子惯运动设置
Sinf.Rlist=[2 4;
            0 0;
             0 0];
Sinf.ulist=[1 1]*arcdeg;
Sinf.flist=[1 1]*0.01;
Sinf.aerrlist=[1 1;
          2 2;
          3 3]*arcdeg;
Sinf.wim0=wim0;
SmoveCell=my_nSmovePack(Sinf);
Smove1=SmoveCell{1};
Smove2=SmoveCell{2};
% Smove1.um=1*arcdeg;%R位置处等效角度振幅
% Smove1.f=0.01;%振动频率
% Smove1.u=[zeros(1,len);Smove1.um*sin(2*pi*Smove1.f*(0:ts:((len-1)*ts)));zeros(1,len)];
% Smove1.R0=[2;0;0];%理想时 主子间0时刻前初始相对位置
% Smove1.Rf0=[0;0;0];%理想时 子0时刻前偏移
% Smove1.Rf=[-Smove1.R0(1)*sin(Smove1.u(2,:)).*sin(Smove1.u(2,:));zeros(1,len);Smove1.R0(1)*sin(Smove1.u(2,:)).*cos(Smove1.u(2,:))];
% Smove1.R=Smove1.R0+Smove1.Rf;
% Smove1.att=2*Smove1.u;
% Smove1.att_s0=[0;0;0]*arcdeg;%理想时 主子间0时刻前初始相对角度
% Smove1.diffRf0=[0;0;0];%Rf的初始变化率
% Smove1.U0=Smove1.diffRf0+cross(wim0,Smove1.R0);
% Smove1.atterr0=[1;2;3]*arcdeg;%注意这里旋转顺序zxy
% %子惯 2 的运动
% Smove2=Smove1;
% Smove2.um=1*arcdeg;%R位置处等效角度振幅
% Smove2.f=0.01;%振动频率
% Smove2.u=[zeros(1,len);Smove2.um*sin(2*pi*Smove2.f*(0:ts:((len-1)*ts)));zeros(1,len)];
% Smove2.R0=[4;0;0];%理想时 主子间0时刻前初始相对位置
% Smove2.Rf0=[0;0;0];%理想时 子0时刻前偏移
% Smove2.Rf=[-Smove2.R0(1)*sin(Smove2.u(2,:)).*sin(Smove2.u(2,:));zeros(1,len);Smove2.R0(1)*sin(Smove2.u(2,:)).*cos(Smove2.u(2,:))];
% Smove2.R=Smove2.R0+Smove2.Rf;
% Smove2.att=2*Smove2.u;
% Smove2.att_s0=[0;0;0]*arcdeg;%理想时 主子间0时刻前初始相对角度
% Smove2.diffRf0=[0;0;0];%Rf的初始变化率
% Smove2.U0=Smove2.diffRf0+cross(wim0,Smove2.R0);
% Smove2.atterr0=[1;2;3]*arcdeg;%注意这里旋转顺序zxy
%% 子惯反解算数据
SINS1=my_invRIpack(MINS,Smove1);
SINS2=my_invRIpack(MINS,Smove2);
%% 角度真值 zxy转序
fprintf('子惯仿真完成！\n纯子惯相对导航...%5.0f %%',0); 
%% 子惯误差设计
SERR1.eb = [0.1;0.1;0.1]*dph; SERR1.web = [0.1;0.1;0.1]*dpsh;   %陀螺常值零偏，角度随机游走
SERR1.db = [20;20;20]*ug; SERR1.wdb = [20;20;20]*ugpsHz;  %加速度计常值偏值，速度随机游走
SERR2.eb = 0*[0.1;0.1;0.1]*dph; SERR2.web = 0*[0.1;0.1;0.1]*dpsh;   %陀螺常值零偏，角度随机游走
SERR2.db = 0*[20;20;20]*ug; SERR2.wdb = 0*[20;20;20]*ugpsHz;  %加速度计常值偏值，速度随机游走
%% 相对导航
atterr0_1=[1.10;2.01;3.50]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
atterr0_2=[1.10;2.01;3.50]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
SINSR1=my_SINSgetResult(MINS,SINS1,Sinf.len/2,SERR1,atterr0_1);
SINSR2=my_SINSgetResult(MINS,SINS2,Sinf.len/2,SERR2,atterr0_2);
fprintf('纯惯相对导航完成！\n滤波...%5.0f %%',0);
%% 滤波
nts = 2*ts; %子样数和采样时间
%子惯1
KFinit1.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit1.rk = [0.0001;0.0001;0.0001];  
KFinit1.Rk = diag(KFinit1.rk)^2;
KFinit1.P0 = diag([[0.1;0.1;10]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR1,~]=my_getFResult(MINS,SINS1,KFinit1,atterr0_1,Sinf.len/2,...
    SINSR1.ws_m_addnoise,SINSR1.fs_m_addnoise);
%子惯2
KFinit2.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit2.rk = [0.0001;0.0001;0.0001];  
KFinit2.Rk = diag(KFinit2.rk)^2;
KFinit2.P0 = diag([[0.1;0.1;10]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR2,~]=my_getFResult(MINS,SINS2,KFinit2,atterr0_2,Sinf.len/2,...
    SINSR2.ws_m_addnoise,SINSR2.fs_m_addnoise);
fprintf('滤波完成！\n'); 
%% 
R=SINS1.R;%理论位置
figure%理论相对位置与纯解算相对位置对比
subplot(311)
T=ts:ts:(Sinf.len*ts);
plot(T(2:4:end),R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.Rall(1:2:end,1),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('x/m');
title('主子惯导相对位置R');
grid on;
subplot(312)
plot(T(2:4:end),R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.Rall(1:2:end,2),'--','LineWidth',2)
grid on;
legend('真实值','计算值');
xlabel('时间/s');ylabel('y/m');
title('主子惯导相对位置R');
subplot(313)
plot(T(2:4:end),R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.Rall(1:2:end,3),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('z/m');
title('主子惯导相对位置R');
grid on;
%%
figure%理论角度与纯相对解算角度对比
subplot(311)
plot(T(2:4:end),(SINS1.atttrue(1,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(1,1:2:end)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('俯仰角/°');
title('主子惯导相对姿态');

grid on;
subplot(312)
plot(T(2:4:end),(SINS1.atttrue(2,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(2,1:2:end)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('滚转角/°');
title('主子惯导相对姿态');

grid on;
subplot(313)
plot(T(2:4:end),(SINS1.atttrue(3,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(3,1:2:end)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('航向角/°');
title('主子惯导相对姿态');

grid on;
%% 
figure%理论相对位置与纯解算相对位置对比
subplot(311)
T=ts:ts:(Sinf.len*ts);
plot(T(2:4:end),SINS1.R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,1),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('x/m');
title('主子惯导相对位置R（滤波）');
grid on;
subplot(312)
plot(T(2:4:end),SINS1.R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,2),'--','LineWidth',2)
grid on;
legend('真实值','计算值');
xlabel('时间/s');ylabel('y/m');
title('主子惯导相对位置R（滤波）');
subplot(313)
plot(T(2:4:end),SINS1.R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,3),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('z/m');
title('主子惯导相对位置R（滤波）');
grid on;
%%
figure%理论角度与纯相对解算角度对比
subplot(311)
plot(T(2:4:end-1),(SINS1.atttrue(1,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(1,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('俯仰角/°');
title('主子惯导相对姿态（滤波）');

grid on;
subplot(312)
plot(T(2:4:end-1),(SINS1.atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(2,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('滚转角/°');
title('主子惯导相对姿态（滤波）');

grid on;
subplot(313)
plot(T(2:4:end-1),(SINS1.atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(3,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('航向角/°');
title('主子惯导相对姿态（滤波）');

grid on;
%% 
figure
subplot(311)
T=ts:ts:(Sinf.len*ts);
plot(T(2:4:end),SINSFR1.Rall(1:2:end,1)'-SINS1.R(1,2:4:end),'LineWidth',2)
legend('X');
xlabel('时间/s');ylabel('x/m');
title('主子惯导相对位置误差');
grid on;
subplot(312)
plot(T(2:4:end),SINSFR1.Rall(1:2:end,2)'-SINS1.R(2,2:4:end),'LineWidth',2)
grid on;
legend('Y');
xlabel('时间/s');ylabel('y/m');
title('主子惯导相对位置误差');
subplot(313)
plot(T(2:4:end),SINSFR1.Rall(1:2:end,3)'-SINS1.R(3,2:4:end),'LineWidth',2)
legend('Z');
xlabel('时间/s');ylabel('z/m');
title('主子惯导相对位置误差');
grid on;
%%
figure
subplot(311)
plot(T(2:4:end-1),(SINSFR1.attall(1,1:2:end-1)-SINS1.atttrue(1,2:4:end-1))/arcdeg,'LineWidth',2)
xlabel('时间/s');ylabel('俯仰角/°');title('主子惯导相对姿态误差');grid on;
subplot(312)
plot(T(2:4:end-1),(SINSFR1.attall(2,1:2:end-1)-SINS1.atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
xlabel('时间/s');ylabel('滚转角/°');title('主子惯导相对姿态误差');grid on;
subplot(313)
plot(T(2:4:end-1),(SINSFR1.attall(3,1:2:end-1)-SINS1.atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
xlabel('时间/s');ylabel('航向角/°');title('主子惯导相对姿态误差');grid on;
%% 统计误差
errlen=1/(Smove1.f*nts);
Rerr=my_getSErr(SINSFR1.Rall(1:end,:)'-SINS1.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSFR1.attall(:,1:end)-SINS1.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSFR1.Rall(1:end-1,:)'-SINS1.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSFR1.attall(:,1:end-1)-SINS1.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('滤波终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1/(Smove1.f*nts);
Rerr=my_getSErr(SINSR1.Rall(1:end,:)'-SINS1.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSR1.attall(:,1:end)-SINS1.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('纯惯统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSR1.Rall(1:end-1,:)'-SINS1.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSR1.attall(:,1:end-1)-SINS1.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('纯惯终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
%%
disp "SINS2"
errlen=1/(Smove2.f*nts);
Rerr=my_getSErr(SINSFR2.Rall(1:end,:)'-SINS2.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSFR2.attall(:,1:end)-SINS2.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSFR2.Rall(1:end-1,:)'-SINS2.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSFR2.attall(:,1:end-1)-SINS2.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('滤波终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1/(Smove2.f*nts);
Rerr=my_getSErr(SINSR2.Rall(1:end,:)'-SINS2.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSR2.attall(:,1:end)-SINS2.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('纯惯统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSR2.Rall(1:end-1,:)'-SINS2.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSR2.attall(:,1:end-1)-SINS2.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('纯惯终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));