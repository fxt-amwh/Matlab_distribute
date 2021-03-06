%% 一主+多子（高度模块化）
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
global setPara;
d2g=180/pi;
g2d=1/d2g;
setPara.ts=0.01;%采样时间
setPara.att0=[0;0;0]*d2g;%载机初始姿态
setPara.vn0=[0;200;0];%载机初始速度
setPara.latitude=34*d2g;%纬度
setPara.longitude=108*d2g;%经度
setPara.h=105;%高
              %俯仰角速率   滚转角速率   航向角速率  加速度     时间
              %  °/s      °/s        °/s        m/s2       s
setPara.wat=[    0,         0,          0,         0,         10       %静止
                 1,         0,          0,         0,         60       %加速
                -1,         0,          0,         0,         60       %匀速
                 0,         0,          0,         0,         50       %匀速
                 0,         0,          1,         0,         360      %匀速转弯
                 0,         0,          0,         0,         50       %匀速
                -1,         0,          0,         0,         60       %匀速
                 1,         0,          0,         0,         60       %减速
                 0,         0,          0,         0,         10       %减速
        ];    %静止

ts = setPara.ts;%采样时间
att0 = setPara.att0; vn0 = setPara.vn0; pos0 = [[setPara.latitude;setPara.longitude];setPara.h];
wat = setPara.wat;
wat(:,1:3) = wat(:,1:3)*arcdeg;  % deg/min->deg/s
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
MINS.wim=wim;MINS.wim0=wim0;MINS.fm=fm;MINS.ts=ts;
fprintf('主惯仿真完成！\n');