%% 一主+多子 尝试将所有子惯数据转换到第一子惯处
% 陀螺的转换只需要旋转
% 加计的转换还需要考虑杆臂加速度
% 版本时间：2021/02/14 一主多子多节点融合初步方案实现
clc
clear
close all
gvar;    % 加载全局变量
ts = 0.05;%采样时间
g2d=180/pi;
d2g=1/g2d;
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;105];
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
MINS.wim=wim;MINS.wim0=wim0;MINS.fm=fm;MINS.ts=ts;
fprintf('主惯仿真完成！\n');
%% 最远子惯信息
simf=0.01;
uf_lw_m=5*d2g;
R_lw=[3;0;0];
%% 全局设计
flag.EnReli=true;%是否纯进行相对导航不滤波
flag.EnFusion=false;%是否进行融合
nameList=["SINS1_1","SINS1_2","SINS1_3","SINS1_4"
    "SINS2_1","SINS2_2","SINS2_3","SINS2_4"];
% 子惯运动数据设计
fprintf('子惯仿真...\n'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;


% Sinf.Rlist=[2 2 2 2 2;
%             0 0 0 0 0;
%             0 0 0 0 0];
Sinf.Rlist=[linspace(1,R_lw(1,1),4);zeros(2,4)];
Sinf.ulist=uf_lw_m*Sinf.Rlist(1,:)/R_lw(1,1);%子惯处挠曲角
% Sinf.Rlist=[R_lw(1,1)*ones(1,4);zeros(2,4)];
% Sinf.ulist=uf_lw_m*ones(1,4);%子惯处挠曲角
Sinf.flist=simf*ones(1,size(Sinf.Rlist,2));%挠曲频率
Sinf.aerrlist=[ 1  1  1  1 ;
                2  2  2  2 ;
                3  3  3  3 ]*d2g;%子惯安装误差角
Sinf.wim0=MINS.wim0;SinfCell{1,1}=Sinf;
% Sinf.Rlist(2,:)=0.5;%第二排y轴位置
% SinfCell{1,2}=Sinf;
Smove=my_nSmovePackN(SinfCell);
[M,N]=size(Smove);%根据设计推断阵面排列
% 子惯误差设计
SERRList.eb=zeros(3,N,M);%zeros(3,N,M)
SERRList.web=zeros(3,N,M);%zeros(3,N,M)
SERRList.db=zeros(3,N,M);%zeros(3,N,M)
SERRList.wdb=zeros(3,N,M);%zeros(3,N,M)
%%
for m=1:M
    for n=1:N
        SERRList.eb(:,n,m)=[1;1;1]/2*dph;%陀螺常值零偏 一列对应一个
        SERRList.web(:,n,m)=[1;1;1]/2*dpsh;%角度随机游走 一列对应一个
        SERRList.db(:,n,m)=[200;200;200]*ug;%加速度计常值偏值 一列对应一个
        SERRList.wdb(:,n,m)=[200;200;200]*ugpsHz;%速度随机游走 一列对应一个
    end
end
% 初始时的安装误差角约值
atterrList=zeros(3,N,M);%zeros(3,N,M);
for m=1:M
    atterrList(:,:,m)=[ 3  1  1  3 ;
                        2  3  2  2 ;
                        3  3  3  3 ]*d2g;%子惯安装误差角
end
%% 子惯反解算数据
SINS=my_invRIpackN(MINS,Smove);
fprintf('子惯仿真完成！\n'); 
%% 子惯误差转换为元组
SERR=my_getSERRN(SERRList);
%% 相对导航
atterr=my_getAtterr(atterrList);
if flag.EnReli
    SINSR=my_SINSgetResultN(MINS,SINS,SinfCell{1,1}.len/2,SERR,atterr); 
end
%% 滤波
flag.sensor=true;%是否量测融合
flag.EnFusion=true;%是否进行融合
nts = 2*ts; %子样数和采样时间
KFinit=my_KFinitN(SERR,SINS,2*ts);
if flag.sensor 
     if flag.EnReli
        % [SINSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,1);%单个子惯分离版
        flag.select=1;
        if flag.EnFusion
            [SINSFR,Filter,NSINSFR]=my_getFResultLoopN_allEn(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,flag);%并行集中版
        else
            [SINSFR,Filter,~]=my_getFResultLoopN_allEn(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,flag);%并行集中版
        end
    else
        flag.select=2;
        % [SINSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,2);%单个子惯分离版
        if flag.EnFusion
            [SINSFR,Filter,NSINSFR]=my_getFResultLoopN_allEn(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,flag);%并行集中版
        else
            [SINSFR,Filter,~]=my_getFResultLoopN_allEn(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,flag);%并行集中版
        end
    end
else
    if flag.EnReli
        % [SINSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,1);%单个子惯分离版
        flag.select=1;
        if flag.EnFusion
            [SINSFR,Filter,NSINSFR]=my_getFResultLoopN_new(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,flag);%并行集中版
        else
            [SINSFR,Filter,~]=my_getFResultLoopN_new(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,flag);%并行集中版
        end
    else
        flag.select=2;
        % [SINSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,2);%单个子惯分离版
        if flag.EnFusion
            [SINSFR,Filter,NSINSFR]=my_getFResultLoopN_new(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,flag);%并行集中版
        else
            [SINSFR,Filter,~]=my_getFResultLoopN_new(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,flag);%并行集中版
        end
    end
end
%% 绘制相对导航子惯1位置结果
flag.Ceil=1;flag.flagSINS_M=1;flag.flagSINS_N=1;flag.figureFlag=1;myfigure;
% 绘制相对导航子惯1姿态结果
flag.figureFlag=2;myfigure;
% 绘制子惯1滤波位置结果
flag.figureFlag=3;myfigure;
% 绘制子惯1滤波姿态结果
flag.figureFlag=4;myfigure;
% 绘制子惯1滤波位置误差
flag.figureFlag=5;myfigure;
% 绘制子惯1滤波姿态误差
flag.figureFlag=6;myfigure;
if flag.EnFusion
    % 绘制子惯1滤波位置误差
    flag.figureFlag=7;myfigure;
%     % 绘制子惯1滤波姿态误差
%     flag.figureFlag=8;myfigure;
end
%% 统计误差
for m=1:M
    for n=1:N
        disp(nameList(m,n));
%         errlen=1/(Smove{m,n}.f*nts);
%         Rerr=my_getSErr(SINSFR{m,n}.Rall(1:end,:)'-SINS{m,n}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
%         Atterr=my_getSErr(SINSFR{m,n}.attall(:,1:end)-SINS{m,n}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
%         fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%             Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
        errlen=1;
        beg=ceil(200/nts);
        End=ceil(500/nts);
        Rerr=my_getSErrBE(SINSFR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),beg,End,2)*1000;        %位置误差单位mm
        Rserr=my_getVErr(SINSFR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),beg,End,2)*1000;        %位置误差单位mm
        Atterr=my_getSErrBE(SINSFR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),beg,End,2)/arcdeg*60; %姿态误差单位分
        Attserr=my_getVErr(SINSFR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),beg,End,2)/arcdeg*60;
        fprintf('滤波终值,位置误差%f %f %f（mm）,位置误差方差%f %f %f(mm)，姿态误差%f %f %f（分）,姿态均方差%f %f %f(分)\n',...
            Rerr(1),Rerr(2),Rerr(3),Rserr(1),Rserr(2),Rserr(3),Atterr(1),Atterr(2),Atterr(3),Attserr(1),Attserr(2),Attserr(3));
%         if flag.EnReli
%             errlen=1/(Smove{m,n}.f*nts);
%             Rerr=my_getSErr(SINSR{m,n}.Rall(1:end,:)'-SINS{m,n}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
%             Atterr=my_getSErr(SINSR{m,n}.attall(:,1:end)-SINS{m,n}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
%             fprintf('纯惯统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%                 Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
%             errlen=1;
%             Rerr=my_getSErr(SINSR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
%             Atterr=my_getSErr(SINSR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
%             fprintf('纯惯终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%                 Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
%         end
    end
end
%%
if flag.EnFusion
    disp('融合结果');
    errlen=1/(Smove{1,1}.f*nts);
    Rerr=my_getSErr(NSINSFR.Rall(1:end,:)'-SINS{1,1}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
    Atterr=my_getSErr(NSINSFR.attnewall(:,1:end)-SINS{1,1}.attuftrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
    fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
        Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
    errlen=1;
    beg=ceil(200/nts);
    End=ceil(500/nts);
    Rerr=my_getSErrBE(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),beg,End,2)*1000;        %位置误差单位mm
    Rserr=my_getVErr(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),beg,End,2)*1000;        %位置误差单位mm
    Atterr=my_getSErrBE(NSINSFR.attnewall(:,1:end-1)-SINS{1,1}.attuftrue(:,1:2:end-2),beg,End,2)/arcdeg*60; %姿态误差单位分
    Attserr=my_getVErr(NSINSFR.attnewall(:,1:end-1)-SINS{1,1}.attuftrue(:,1:2:end-2),beg,End,2)/arcdeg*60;
    fprintf('滤波终值,位置误差%f %f %f（mm）,位置误差方差%f %f %f(mm)，姿态误差%f %f %f（分）,姿态均方差%f %f %f(分)\n',...
        Rerr(1),Rerr(2),Rerr(3),Rserr(1),Rserr(2),Rserr(3),Atterr(1),Atterr(2),Atterr(3),Attserr(1),Attserr(2),Attserr(3));   
end
%% 绘图
atttrueFR=NSINSFR.attnewall(:,1:end-1);
atttrue=SINS{1,1}.attuftrue(:,1:2:end-2);
T=ts:ts:(SinfCell{1,1}.len*ts);
figure
subplot(311)
plot(T(2:2:end-1),atttrueFR(1,:)-atttrue(1,:),'LineWidth',2)
xlabel('时间/s');ylabel('俯仰角/°');title('主子惯导相对姿态误差（融合）');grid on;ylim([-0.01,0.01]);
subplot(312)
plot(T(2:2:end-1),atttrueFR(2,:)-atttrue(2,:),'LineWidth',2)
xlabel('时间/s');ylabel('滚转角/°');title('主子惯导相对姿态误差（融合）');grid on;ylim([-0.01,0.01]);
subplot(313)
plot(T(2:2:end-1),atttrueFR(3,:)-atttrue(3,:),'LineWidth',2)
xlabel('时间/s');ylabel('航向角/°');title('主子惯导相对姿态误差（融合）');grid on;ylim([-0.01,0.01]);

figure%理论角度与纯相对解算角度对比
subplot(311)
plot(T(2:2:end-1),atttrue(1,:),'LineWidth',2)
hold on;
plot(T(2:2:end-1),atttrueFR(1,:),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('俯仰角/°');
title('主子惯导相对姿态（融合）');
grid on

subplot(312)
plot(T(2:2:end-1),atttrue(2,:),'LineWidth',2)
hold on;
plot(T(2:2:end-1),atttrueFR(2,:),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('滚转角/°');
title('主子惯导相对姿态（融合）');
grid on;

subplot(313)
plot(T(2:2:end-1),atttrue(3,:),'LineWidth',2)
hold on;
plot(T(2:2:end-1),atttrueFR(3,:),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('航向角/°');
title('主子惯导相对姿态（融合）');
grid on