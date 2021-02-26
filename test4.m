%% 一主+3子（中间备份）
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
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;105];
%     俯仰角速率 横滚角速率 方位角速率 纵向加速度 持续时间 （环形轨迹设计）
wat = [  0,         -60,          120,         0,         30               %30
        0,         60,          -120,        0,         60                 %90
        0,         0,          0,         0,         120                   %210
        -60,         0,          0,         0,         30                  %240
%         0,         0,          360*60,         0,         360       %匀速转弯
        0,         0,          0,         0,         120                   %360
        60,         0,          0,         0,         30                   %390
        0,         0,          0,         0,         220                   %610
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
MINS.ts=ts;
fprintf('主惯仿真完成！\n');
%% 子惯运动数据设计
fprintf('子惯仿真...'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;
Sinf.Rlist=[15 15 15;
            0 0 0;
            0 0 0];
Sinf.ulist=[1 2 3]*arcdeg;
Sinf.flist=[1 1 1]*0.01;
Sinf.aerrlist=[ 1 1 1;
                2 2 2;
                3 3 3]*arcdeg;
Sinf.wim0=MINS.wim0;
SmoveCell=my_nSmovePack(Sinf);
Smove1=SmoveCell{1};
Smove2=SmoveCell{2};
Smove3=SmoveCell{3};
%% 子惯反解算数据
SINS1=my_invRIpack(MINS,Smove1);
SINS2=my_invRIpack(MINS,Smove2);
SINS3=my_invRIpack(MINS,Smove3);
%% 角度真值 zxy转序
fprintf('子惯仿真完成！\n纯子惯相对导航...%5.0f %%',0); 
%% 子惯误差设计
SERR1.eb = [0.1;0.1;0.1]*dph; SERR1.web = [0.1;0.1;0.1]*dpsh;   %陀螺常值零偏，角度随机游走
SERR1.db = [20;20;20]*ug; SERR1.wdb = [20;20;20]*ugpsHz;  %加速度计常值偏值，速度随机游走
SERR2.eb = [0.1;0.1;0.1]*dph; SERR2.web = [0.1;0.1;0.1]*dpsh;   %陀螺常值零偏，角度随机游走
SERR2.db = [20;20;20]*ug; SERR2.wdb = [20;20;20]*ugpsHz;  %加速度计常值偏值，速度随机游走
SERR3.eb = [0.1;0.1;0.1]*dph; SERR3.web = [0.1;0.1;0.1]*dpsh;   %陀螺常值零偏，角度随机游走
SERR3.db = [20;20;20]*ug; SERR3.wdb = [20;20;20]*ugpsHz;  %加速度计常值偏值，速度随机游走
%% 相对导航
atterr0_1=[1.00;2.00;3.00]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
atterr0_2=[1.00;2.00;3.00]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
atterr0_3=[1.00;2.00;3.00]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
SINSR1=my_SINSgetResult(MINS,SINS1,Sinf.len/2,SERR1,atterr0_1);
SINSR2=my_SINSgetResult(MINS,SINS2,Sinf.len/2,SERR2,atterr0_2);
SINSR3=my_SINSgetResult(MINS,SINS3,Sinf.len/2,SERR3,atterr0_3);
fprintf('\n纯惯相对导航完成！\n滤波...%5.0f %%',0);  
%% 滤波
nts = 2*ts; %子样数和采样时间
%子惯1
KFinit1.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit1.rk = [0.001;0.001;0.001]*sqrt(SINS1.R0'*SINS1.R0);  
KFinit1.Rk = diag(KFinit1.rk)^2;
KFinit1.P0 = diag([[0.1;0.1;0.1]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR1,Filter1]=my_getFResult(MINS,SINS1,KFinit1,atterr0_1,Sinf.len/2,...
    SINSR1.ws_m_addnoise,SINSR1.fs_m_addnoise);
%子惯2
KFinit2.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit2.rk = [0.001;0.001;0.001]*sqrt(SINS2.R0'*SINS2.R0);   
KFinit2.Rk = diag(KFinit2.rk)^2;
KFinit2.P0 = diag([[0.1;0.1;1]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR2,Filter2]=my_getFResult(MINS,SINS2,KFinit2,atterr0_2,Sinf.len/2,...
    SINSR2.ws_m_addnoise,SINSR2.fs_m_addnoise);
%子惯3
KFinit3.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit3.rk = [0.001;0.001;0.001]*sqrt(SINS3.R0'*SINS3.R0);  
KFinit3.Rk = diag(KFinit3.rk)^2;
KFinit3.P0 = diag([[0.1;0.1;10]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR3,Filter3]=my_getFResult(MINS,SINS3,KFinit3,atterr0_3,Sinf.len/2,...
    SINSR3.ws_m_addnoise,SINSR3.fs_m_addnoise);
fprintf('\n滤波完成！\n'); 
%% 绘制相对导航子惯1位置结果
flag.Ceil=0;flag.flagSINS_M=1;flag.flagSINS_N=1;flag.figureFlag=1;
myfigure;
%% 绘制相对导航子惯1姿态结果
flag.figureFlag=2;myfigure;
%% 绘制子惯1滤波位置结果
flag.figureFlag=3;myfigure;
%% 绘制子惯1滤波姿态结果
flag.figureFlag=4;myfigure;
%% 绘制子惯1滤波位置误差
flag.figureFlag=5;myfigure;
%% 绘制子惯1滤波姿态误差
flag.figureFlag=6;myfigure;
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
%%
disp "SINS3"
errlen=1/(Smove3.f*nts);
Rerr=my_getSErr(SINSFR3.Rall(1:end,:)'-SINS3.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSFR3.attall(:,1:end)-SINS3.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSFR3.Rall(1:end-1,:)'-SINS3.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSFR3.attall(:,1:end-1)-SINS3.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('滤波终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1/(Smove2.f*nts);
Rerr=my_getSErr(SINSR3.Rall(1:end,:)'-SINS3.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSR3.attall(:,1:end)-SINS3.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('纯惯统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSR3.Rall(1:end-1,:)'-SINS3.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
Atterr=my_getSErr(SINSR3.attall(:,1:end-1)-SINS3.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
fprintf('纯惯终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));