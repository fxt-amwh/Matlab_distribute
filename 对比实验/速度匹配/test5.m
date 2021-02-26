%文件名称：test5
%作者信息：冯鑫涛
%功能描述：获得主惯数据结构体以及子惯数据结构体 不保存 单纯为了方便测试跟踪内部数据
%版本时间：2021/1/30 19:17
%% 主惯运动数据设计仿真
clc
clear
close all
gvar;    % 加载全局变量
ts = 0.01;%采样时间
g2d=180/pi;
d2g=1/g2d;
att0 = [0;0;0]*arcdeg; vn0 = [0;100;0]; pos0 = [[34;108]*arcdeg;100];
%     俯仰角速率 横滚角速率 方位角速率 纵向加速度 持续时间 （环形轨迹设计）
% wat = [  0,         -60,          120,         0,         30               %30
%         0,         60,          -120,        0,         60                 %90
%         0,         0,          0,         0,         120                   %210
%         -60,         0,          0,         0,         30                  %240
%         0,         0,          0,         0,         120                   %360
%         60,         0,          0,         0,         30                   %390
%         0,         0,          0,         0,         220                   %610
%         ];    %静止
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
% 组合TR
TR.att0=att0;
TR.vn0=vn0;
TR.pos0=pos0;
TR.vn=vn;
TR.attm=attm;
TR.pos=pos;
%% 主惯信息
wim0=my_getII(wm(1,:)',ts);
wim=[wim0,my_getII(wm',ts)];%wim注意是在主坐标系下的量测
fm=[my_getII(vm(1,:)',ts),my_getII(vm',ts)];
MINS.wim=wim;MINS.wim0=wim0;MINS.fm=fm;MINS.ts=ts;
fprintf('主惯仿真完成！\n');
%% 最远子惯信息
simf=0.01;
uf_lw_m=5*d2g;
R_lw=[2;0;0];
%% 全局设计
flag.EnReli=false;%是否纯进行相对导航不滤波
flag.EnFusion=true;%是否进行融合
nameList=["SINS1_1","SINS1_2"];
[M,N]=size(nameList);%根据设计推断阵面排列
% 子惯运动数据设计
fprintf('子惯仿真...\n'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;


Sinf.Rlist=[R_lw(1,1)*ones(1,N);zeros(2,N)];
Sinf.ulist=uf_lw_m*Sinf.Rlist(1,:)/R_lw(1,1);%子惯处挠曲角
Sinf.ulist(:,:)=0;%子惯处挠曲角
Sinf.flist=simf*ones(1,size(Sinf.Rlist,2));%挠曲频率
Sinf.aerrlist=0*[ 1  1 ;
                2  2 ;
                3  3 ]*d2g;%子惯安装误差角
Sinf.wim0=MINS.wim0;SinfCell{1,1}=Sinf;
 
Smove=my_nSmovePackN(SinfCell);
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
%         SERRList.eb(:,n,m)=[0.1;0.1;0.1]*dph;%陀螺常值零偏 一列对应一个
%         SERRList.web(:,n,m)=[0.1;0.1;0.1]*dpsh;%角度随机游走 一列对应一个
%         SERRList.db(:,n,m)=[20;20;20]*ug;%加速度计常值偏值 一列对应一个
%         SERRList.wdb(:,n,m)=[20;20;20]*ugpsHz;%速度随机游走 一列对应一个
%         SERRList.eb(:,n,m)=[0.5;0.5;0.5]*dph;%陀螺常值零偏 一列对应一个
%         SERRList.web(:,n,m)=[0.5;0.5;0.5]*dpsh;%角度随机游走 一列对应一个
%         SERRList.db(:,n,m)=[100;100;100]*ug;%加速度计常值偏值 一列对应一个
%         SERRList.wdb(:,n,m)=[100;100;100]*ugpsHz;%速度随机游走 一列对应一个
        SERRList.eb(:,n,m)=[1;1;1]*dph;%陀螺常值零偏 一列对应一个
        SERRList.web(:,n,m)=[1;1;1]*dpsh;%角度随机游走 一列对应一个
        SERRList.db(:,n,m)=[200;200;200]*ug;%加速度计常值偏值 一列对应一个
        SERRList.wdb(:,n,m)=[200;200;200]*ugpsHz;%速度随机游走 一列对应一个
    end
end
% 初始时的安装误差角约值
atterrList=zeros(3,N,M);%zeros(3,N,M);
for m=1:M
    atterrList(:,:,m)=0*[ 1  1 ;
                        2  2 ;
                        3  3 ]*d2g;%子惯安装误差角
end
%% 子惯反解算数据
SINS=my_invRIpackN(MINS,Smove);
fprintf('子惯仿真完成！\n'); 
%% 子惯误差转换为元组
SERR=my_getSERRN(SERRList);
%% 子惯误差注入
len=length(SINS{1,2}.wis);
for i=1:len/2
    [ws_m, fs_m] = my_imuadderr(SINS{1,2}.wis(:,(2*i-1):(2*i))', SINS{1,2}.fs(:,(2*i-1):(2*i))', SERR{1,2}.eb, SERR{1,2}.web, SERR{1,2}.db, SERR{1,2}.wdb, ts);%子惯注入噪声
end
%%
% qnb=a2qua(att0);
% vn=vn0;
% posSINS=pos0;
% % posSINS(3)=110;
% len=length(vm)-1;
% nn=2;nts=nn*ts;
% SINS=zeros(floor(len/2),10);kk=1;t=0;
% % eb=[0.01;0.01;0.01]*dph;
% % web=[0.001;0.001;0.001]*dpsh;
% % db=[10;10;10]*ug;
% % wdb=[0;0;0]*ugpsHz;
% eb=[0;0;0]*dph;
% web=[0;0;0]*dpsh;
% db=[0;0;0]*ug;
% wdb=[0;0;0]*ugpsHz;
% sensor1=zeros(2,6,len/2);
% sensor=zeros(2,6,len/2);
% num=1;
% for k=1:nn:len
%     t=t+nts;
%     [wm1,vm1]=imuadderr(wm(k:k+1,:),vm(k:k+1,:), eb, web, db, wdb, ts);
%     sensor1(:,:,num)=[wm1,vm1];sensor(:,:,num)=[wm1,vm1];
%     num=num+1;
%     [qnb,vn,posSINS]=my_insupdate(qnb,vn,posSINS,wm1,vm1,ts);
%     SINS(kk,:)=[q2att(qnb);vn;posSINS;t];kk=kk+1;
%     
%     if mod(t,100)<nts,disp(fix(t));end
% end
% %
% figure
% delpos=deltapos(pos);
% plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
% hold on;
% delposSINS=deltapos(SINS(:,7:9));
% plot3(delposSINS(:,1),delposSINS(:,2),delposSINS(:,3));
% grid on;
% xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');
% err=delposSINS(end,:)-delpos(end,:);
% disp '纬度误差 经度误差 高度误差'
% disp(err)