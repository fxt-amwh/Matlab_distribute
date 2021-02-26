%% 一主+多子（高度模块化）使用角度分别修正后的多节点融合
% 1、根据设计轨迹运动反解算主惯传感器输出
% 2、根据设计机翼运动反解算子惯传感器数据
% 3、子惯数据注入误差
% 4、纯惯性相对导航
% 5、滤波
% 6、误差计算
% 版本时间：2020/12/19 多节点初步测试
%% 主惯运动数据设计仿真
clc
clear
close all
gvar;    % 加载全局变量
ts = 0.01;%采样时间
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;105];
%     俯仰角速率 横滚角速率 方位角速率 纵向加速度 持续时间 （环形轨迹设计）
% wat = [  0,         -60,          120,         0,         30               %30
%         0,         60,          -120,        0,         60                 %90
%         0,         0,          0,         0,         120                   %210
%         -60,         0,          0,         0,         30                  %240
%         0,         0,          0,         0,         120                   %360
%         60,         0,          0,         0,         30                   %390
%         0,         0,          0,         0,         220                   %610
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
%% 全局设计
flag.EnReli=false;%是否纯进行相对导航不滤波
flag.EnFusion=true;%是否进行融合
nameList=["SINS1","SINS2","SINS3","SINS4"];
% 子惯运动数据设计
fprintf('子惯仿真...'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;
Sinf.Rlist=[1 2 3 4;
            0 0 0 0;
            0 0 0 0];
Sinf.ulist=[1 2 3 4]*arcdeg;%子惯处挠曲角
Sinf.flist=[1 1 1 1]*0.01;%挠曲频率
% Sinf.aerrlist=[ 1 1 1 1;
%                 2 2 2 2;
%                 3 3 3 3]*arcdeg;%子惯安装误差角
Sinf.aerrlist=[ 0 0 0 0;
                0 0 0 0;
                0 0 0 0]*arcdeg;%子惯安装误差角
Sinf.wim0=MINS.wim0;SinfCell{1,1}=Sinf;
Smove=my_nSmovePackN(SinfCell);
[M,N]=size(Smove);%根据设计推断阵面排列
% 子惯误差设计
SERRList.eb=zeros(3,N,M);%zeros(3,N,M)
SERRList.web=zeros(3,N,M);%zeros(3,N,M)
SERRList.db=zeros(3,N,M);%zeros(3,N,M)
SERRList.wdb=zeros(3,N,M);%zeros(3,N,M)
% for m=1:M
%     SERRList.eb(:,:,m)=[0.1 0.1 0.1;%
%                         0.1 0.1 0.1;%
%                         0.1 0.1 0.1]*dph;%陀螺常值零偏 一列对应一个
%     SERRList.web(:,:,m)=[0.1 0.1 0.1;%
%                          0.1 0.1 0.1;%
%                          0.1 0.1 0.1]*dpsh;%角度随机游走 一列对应一个
%     SERRList.db(:,:,m)=[20 20 20;%
%                         20 20 20;%
%                         20 20 20]*ug;%加速度计常值偏值 一列对应一个
%     SERRList.wdb(:,:,m)=[20 20 20;%
%                          20 20 20;%
%                          20 20 20]*ugpsHz;%速度随机游走 一列对应一个
% end
for m=1:M
    for n=1:N
        SERRList.eb(:,n,m)=[1;1;1]*dph;%陀螺常值零偏 一列对应一个
        SERRList.web(:,n,m)=[1;1;1]*dpsh;%角度随机游走 一列对应一个
        SERRList.db(:,n,m)=[200;200;200]*ug;%加速度计常值偏值 一列对应一个
        SERRList.wdb(:,n,m)=[200;200;200]*ugpsHz;%速度随机游走 一列对应一个
    end
end
% 初始时的安装误差角约值
atterrList=zeros(3,N,M);%zeros(3,N,M);
for m=1:M
%     atterrList(:,:,m)=[ 1.01 1.01 1.01 1.01;%注意这里旋转顺序zxy 不准确的初始角
%                        2.01 2.01 2.01 2.01;%注意这里旋转顺序zxy 不准确的初始角
%                        3.01 3.01 3.01 3.01]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
    atterrList(:,:,m)=[ 0.01 0.01 0.01 0.01;%注意这里旋转顺序zxy 不准确的初始角
                        0.01 0.01 0.01 0.01;%注意这里旋转顺序zxy 不准确的初始角
                        0.01 0.01 0.01 0.01]*arcdeg;%注意这里旋转顺序zxy 不准确的初始角
%     for n=1:N
%         Smove{m,n}.atterrList=[0;0;0];
%     end
                   
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
nts = 2*ts; %子样数和采样时间
KFinit=my_KFinitN(SERR,SINS,2*ts);
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
    % 绘制子惯1滤波姿态误差
    flag.figureFlag=8;myfigure;
end

%% 统计误差
for m=1:M
    for n=1:N
        disp(nameList(m,n));
        errlen=1/(Smove{m,n}.f*nts);
        Rerr=my_getSErr(SINSFR{m,n}.Rall(1:end,:)'-SINS{m,n}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
        Atterr=my_getSErr(SINSFR{m,n}.attall(:,1:end)-SINS{m,n}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
        fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
            Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
        errlen=1;
        beg=ceil(200/nts);
        End=ceil(500/nts);
        Rerr=my_getSErr(SINSFR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
        Rserr=my_getVErr(SINSFR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),beg,End,2)*1000;        %位置误差单位mm
        Atterr=my_getSErr(SINSFR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
        Attserr=my_getVErr(SINSFR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),beg,End,2)/arcdeg*60;
        fprintf('滤波终值,位置误差%f %f %f（mm）,位置误差方差%f %f %f(mm)，姿态误差%f %f %f（分）,姿态均方差%f %f %f(分)\n',...
            Rerr(1),Rerr(2),Rerr(3),Rserr(1),Rserr(2),Rserr(3),Atterr(1),Atterr(2),Atterr(3),Attserr(1),Attserr(2),Attserr(3));
        if flag.EnReli
            errlen=1/(Smove{m,n}.f*nts);
            Rerr=my_getSErr(SINSR{m,n}.Rall(1:end,:)'-SINS{m,n}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
            Atterr=my_getSErr(SINSR{m,n}.attall(:,1:end)-SINS{m,n}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
            fprintf('纯惯统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
                Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
            errlen=1;
            Rerr=my_getSErr(SINSR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
            Atterr=my_getSErr(SINSR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
            fprintf('纯惯终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
                Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
        end
    end
end
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
    Rerr=my_getSErr(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
    Rserr=my_getVErr(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),beg,End,2)*1000;        %位置误差单位mm
    Atterr=my_getSErr(NSINSFR.attnewall(:,1:end-1)-SINS{1,1}.attuftrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
    Attserr=my_getVErr(NSINSFR.attnewall(:,1:end-1)-SINS{1,1}.attuftrue(:,1:2:end-2),beg,End,2)/arcdeg*60;
    fprintf('滤波终值,位置误差%f %f %f（mm）,位置误差方差%f %f %f(mm)，姿态误差%f %f %f（分）,姿态均方差%f %f %f(分)\n',...
        Rerr(1),Rerr(2),Rerr(3),Rserr(1),Rserr(2),Rserr(3),Atterr(1),Atterr(2),Atterr(3),Attserr(1),Attserr(2),Attserr(3));   
end