%% 一主 + 一子 重要
%文件名称：test11
%作者信息：冯鑫涛
%功能描述：相对导航滤波 有反馈
%版本时间：2021/1/30 19:15
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
figure
delpos=deltapos(pos);
plot3(delpos(:,1),delpos(:,2),delpos(:,3));
% axis([0,-8000,]);
grid on;
xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');
%% 主惯信息
wim0=my_getII(wm(1,:)',ts);
wim=[wim0,my_getII(wm',ts)];%wim注意是在主坐标系下的量测
fm=[my_getII(vm(1,:)',ts),my_getII(vm',ts)];
MINS.wim=wim;
MINS.wim0=wim0;
MINS.fm=fm;
fprintf('主惯仿真完成！\n');
%% 子惯1运动数据设计
fprintf('子惯仿真...'); 
len=size(wm,1)+1;% len=7200;
um=10*arcdeg;
f=0.01;
u=[zeros(1,len);um*sin(2*pi*f*(0:ts:((len-1)*ts)));zeros(1,len)];
R0=[2;0;0];%理想时 主子间0时刻前初始相对位置
Rf0=[0;0;0];%理想时 子0时刻前偏移
Rf=[-R0(1)*sin(u(2,:)).*sin(u(2,:));zeros(1,len);R0(1)*sin(u(2,:)).*cos(u(2,:))];
R=R0+Rf;
att=2*u;
att_s0=[0;0;0]*arcdeg;%理想时 主子间0时刻前初始相对角度
%% 子惯1反解算数据
% atterr0=[1;2;3]*arcdeg;%注意这里旋转顺序zxy
atterr0=[0;0;0]*arcdeg;%注意这里旋转顺序zxy
Cmserr=my_a2mat(atterr0);%注意中间处理旋转顺序yxz
diffRf0=[0;0;0];%Rf的初始变化率
U0=diffRf0+cross(wim0,R0);
SINS1=my_invRI(R,att,atterr0,wim,fm,ts,R0,att_s0,U0);
fprintf('子惯仿真完成！\n纯子惯相对导航...%5.0f %%',0);  
%% 子惯2运动数据设计
fprintf('子惯仿真...'); 
len=size(wm,1)+1;% len=7200;
um=20*arcdeg;
u=[zeros(1,len);um*sin(2*pi*f*(0:ts:((len-1)*ts)));zeros(1,len)];
R0=[4;0;0];%理想时 主子间0时刻前初始相对位置
Rf0=[0;0;0];%理想时 子0时刻前偏移
Rf=[-R0(1)*sin(u(2,:)).*sin(u(2,:));zeros(1,len);R0(1)*sin(u(2,:)).*cos(u(2,:))];
R2=R0+Rf;
att2=2*u;
att_s0=[0;0;0]*arcdeg;%理想时 主子间0时刻前初始相对角度
%% 子惯2反解算数据
% atterr0=[1;2;3]*arcdeg;%注意这里旋转顺序zxy
atterr0=[0;0;0]*arcdeg;%注意这里旋转顺序zxy
Cmserr=my_a2mat(atterr0);%注意中间处理旋转顺序yxz
diffRf0=[0;0;0];%Rf的初始变化率
U0=diffRf0+cross(wim0,R0);
SINS2=my_invRI(R2,att2,atterr0,wim,fm,ts,R0,att_s0,U0);
fprintf('子惯仿真完成！\n纯子惯相对导航...%5.0f %%',0);  
%% 子惯误差设计
MERR.eb=[1;1;1]*dph;
MERR.web=[1;1;1]*dpsh;
MERR.db=[200;200;200]*ug;
MERR.wdb=[200;200;200]*ugpsHz;
SERR1.eb = [1;1;1]*dph; SERR1.web = [1;1;1]*dpsh;   %陀螺常值零偏，角度随机游走
SERR1.db = [200;200;200]*ug; SERR1.wdb = [200;200;200]*ugpsHz;  %加速度计常值偏值，速度随机游走
% SERR1.eb = 0*[1;1;1]*dph; SERR1.web = 0*[1;1;1]*dpsh;   %陀螺常值零偏，角度随机游走
% SERR1.db = 0*[200;200;200]*ug; SERR1.wdb = 0*[200;200;200]*ugpsHz;  %加速度计常值偏值，速度随机游走
%% 相对导航
atterr0_1=[1;0;0]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
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
%% 
R=SINS1.R;%理论位置
figure%理论相对位置与纯解算相对位置对比
subplot(311)
T=ts:ts:(len*ts);
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
plot(T(2:4:end),(atttrue(1,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(1,1:2:end)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('俯仰角/°');
title('主子惯导相对姿态');

grid on;
subplot(312)
plot(T(2:4:end),(atttrue(2,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(2,1:2:end)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('滚转角/°');
title('主子惯导相对姿态');

grid on;
subplot(313)
plot(T(2:4:end),(atttrue(3,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(3,1:2:end)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('航向角/°');
title('主子惯导相对姿态');

grid on;
%% 
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
xlabel('时间/s');ylabel('俯仰角/°');title('主子惯导相对姿态误差');grid on;
subplot(312)
plot(T(2:4:end-1),(SINSFR1.attall(2,1:2:end-1)-atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
xlabel('时间/s');ylabel('滚转角/°');title('主子惯导相对姿态误差');grid on;
subplot(313)
plot(T(2:4:end-1),(SINSFR1.attall(3,1:2:end-1)-atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
xlabel('时间/s');ylabel('航向角/°');title('主子惯导相对姿态误差');grid on;
% %% 统计误差
% errlen=1/(f*nts);
% Rerr=my_getSErr(SINSFR1.Rall(1:end,:)'-R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
% Atterr=my_getSErr(SINSFR1.attall(:,1:end)-atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
% fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%     Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
% errlen=1;
% Rerr=my_getSErr(SINSFR1.Rall(1:end-1,:)'-R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
% Atterr=my_getSErr(SINSFR1.attall(:,1:end-1)-atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
% fprintf('滤波终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%     Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
% errlen=1/(f*nts);
% Rerr=my_getSErr(SINSR1.Rall(1:end,:)'-R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
% Atterr=my_getSErr(SINSR1.attall(:,1:end)-atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
% fprintf('纯惯统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%     Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
% errlen=1;
% Rerr=my_getSErr(SINSR1.Rall(1:end-1,:)'-R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
% Atterr=my_getSErr(SINSR1.attall(:,1:end-1)-atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
% fprintf('纯惯终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%     Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
%% 统计MSE
beg_T=200;
end_T=500;
beg_index=beg_T/nts;
end_index=end_T/nts;
%diag(resultR(:,sam.begin:sam.end)*resultR(:,sam.begin:sam.end)')/length(resultR(:,sam.begin:sam.end))/simLen^2
errall=SINSFR1.attall(:,1:end)-atttrue(:,1:2:end);
err=errall(:,beg_index:end_index);
MSE = 1000*sqrt(diag(err*err'/length(err)));
fprintf("MSE:(mrad)");
disp(MSE');