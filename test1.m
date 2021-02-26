%%
%功能：相对导航解算主惯+子惯
%% 主惯数据仿真
clc
clear
close all
gvar;    % 加载全局变量
ts = 0.1;
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;105];
%     俯仰角速率 横滚角速率 方位角速率 纵向加速度 持续时间 （环形轨迹设计）
wat = [  0,         0,          0,         0,         10       %静止
        60,         0,          0,         0,         60       %加速
        -60,         0,          0,         0,         60       %加速
        0,         0,          0,         0,         50        %匀速
        0,         0,          60,         0,         360       %匀速转弯
        0,         0,          0,         0,         50        %匀速
        -60,         0,          0,         0,         60       %加速
        60,         0,          0,         0,         60        %减速
        0,         0,          0,         0,         10        %减速
        ];    %静止
wat(:,1:3) = wat(:,1:3)*arcdeg/60;  % ->deg/s
[attm, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
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
plot3(delpos(:,1),delpos(:,2),delpos(:,3));
% axis([0,-8000,]);
grid on;
xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');
%% 子惯运动数据
len=size(wm,1)+1;
% len=7200;
um=0.05;
f=0.01;
u=[zeros(1,len);um*sin(2*pi*f*(0:ts:((len-1)*ts)));zeros(1,len)];
R0=[5;0;0];
Rf0=[0;0;0];
Rf=[-R0(1)*sin(u(2,:)).*sin(u(2,:));zeros(1,len);R0(1)*sin(u(2,:)).*cos(u(2,:))];
R=R0+Rf;
att=2*u;
att0=[0;0;0]*arcdeg;
%% 主惯信息
wim0=my_getII(wm(1,:)',ts);
wim=[wim0,my_getII(wm',ts)];%wim注意是在主坐标系下的量测
fm=[my_getII(vm(1,:)',ts),my_getII(vm',ts)];
% wim=zeros(3,len);%wim注意是在主坐标系下的量测
% wim0=[0;0;0];
% for i=1:len
%     wim(:,i)=[6.04543744001202e-05;1.05879118406790e-20;4.07769904065001e-05]*2000;
% end
% wim0=[6.04543744001202e-05;1.05879118406790e-20;4.07769904065001e-05]*2000;
% fm=zeros(3,len);
%% 子惯反解算数据
atterr0=[1;2;3]*arcdeg;%注意这里旋转顺序zxy
diffRf0=[0;0;0];
U0=diffRf0+cross(wim0,R0);
% [wis,fs,U,Cms_cell]=my_invRI(R,att,atterr0,wim,fm,ts,R0,att0,U0);
% SINS1=my_invRI2(R,att,attm',atterr0,wim,fm,ts,R0,att0,U0);
SINS1=my_invRI(R,att,atterr0,wim,fm,ts,R0,att0,U0);

eb = [0.1;0.1;0.1]*dph; web = [0.1;0.1;0.1]*dpsh;   %陀螺常值零偏，角度随机游走
db = [20;20;20]*ug; wdb = [20;20;20]*ugpsHz;  %加速度计常值偏值，速度随机游走

%% 相对导航
% atterr0=[1.01;2;3]*deg2arc;%注意这里旋转顺序zxy
Cmserr=my_a2mat(atterr0);%注意中间处理旋转顺序yxz
Rloop=R0;
Uloop=cross(wim0,R0);
qms=m2qua(Cmserr);

Rall=zeros(len/2,3);
vnall=zeros(len/2,3);
qmsall=zeros(len/2,4);
attall=zeros(3,len/2);

for i=1:(len/2-1)
    [ws_m, fs_m] = my_imuadderr(SINS1.wis(:,(2*i-1):(2*i))', SINS1.fs(:,(2*i-1):(2*i))', eb, web, db, wdb, ts);%子惯注入噪声
    wm_m=wim(:,(2*i-1):(2*i))'; fm_m=fm(:,(2*i-1):(2*i))';

    [qms,vn,Uloop,Rloop] = my_relinsupdate5(qms,Uloop,Rloop,wm_m, fm_m, ws_m, fs_m, ts);% 
%     wm_m=wim(:,(2*i-1):(2*i+1))'; fm_m=fm(:,(2*i-1):(2*i+1))';
%     ws_m=SINS1.wis(:,(2*i-1):(2*i+1))';fs_m=SINS1.fs(:,(2*i-1):(2*i+1))';
% 
%     [qms,vn,Uloop,Rloop] = my_relinsupdate6(qms,Uloop,Rloop,wm_m, fm_m, ws_m, fs_m, ts);
    qmsall(i,:)=qms';
    Rall(i,:)=Rloop';
    vnall(i,:)=vn';
end
for i=1:length(qmsall)
%     attall(:,i)=q2att(qmsall(i,:));
    [rz,rx,ry]=dcm2angle(q2mat(qmsall(i,:)),'zxy');%转换到传统的坐标
    attall(:,i)=[rx,ry,-rz];
end
% [rz,rx,ry]=dcm2angle(a2mat(atterr0),'zxy');
% atterr0=-[rx,ry,-rz]';
%% 
figure
subplot(311)
T=ts:ts:(len*ts);
plot(T(2:4:end),R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),Rall(1:2:end,1),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('x/m');
title('主子惯导相对位置R');
axis([0,800,R0(1)-0.1,R0(1)+0.002]);
grid on;
subplot(312)
plot(T(2:4:end),R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),Rall(1:2:end,2),'--','LineWidth',2)
axis([0,800,-0.003,0.003]);
grid on;
legend('真实值','计算值');
xlabel('时间/s');ylabel('y/m');
title('主子惯导相对位置R');
subplot(313)
plot(T(2:4:end),R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),Rall(1:2:end,3),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('z/m');
title('主子惯导相对位置R');
grid on;
%%
atttrue=att;
for i=1:len
%    atttrue(:,i)=m2att(Cms_cell{i}); 
   [rz,rx,ry]=dcm2angle(SINS1.Cms{i},'zxy');
   atttrue(:,i)=[rx,ry,-rz]; 
end
%%
figure
subplot(311)
plot(T(2:4:end-1),(atttrue(1,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),attall(1,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('俯仰角/°');
title('主子惯导相对姿态');
axis([0,800,atterr0(1)/arcdeg-0.2*um*sqrt(R0'*R0)/arcdeg,atterr0(1)/arcdeg+0.2*um*sqrt(R0'*R0)/arcdeg]);
grid on;
subplot(312)
plot(T(2:4:end-1),(atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),attall(2,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('滚转角/°');
title('主子惯导相对姿态');
axis([0,800,atterr0(2)/arcdeg-um*sqrt(R0'*R0)/arcdeg,atterr0(2)/arcdeg+um*sqrt(R0'*R0)/arcdeg]);
grid on;
subplot(313)
plot(T(2:4:end-1),(atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),attall(3,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('航向角/°');
title('主子惯导相对姿态');
axis([0,800,atterr0(3)/arcdeg-0.05*um*sqrt(R0'*R0)/arcdeg,atterr0(3)/arcdeg+0.05*um*sqrt(R0'*R0)/arcdeg]);
grid on;