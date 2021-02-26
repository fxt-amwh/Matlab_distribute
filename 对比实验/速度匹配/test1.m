%文件名称：test1
%作者信息：冯鑫涛
%功能描述：主子惯导分别进行导航得到速度
%版本时间：2021/1/27 17:19
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
%% 主惯纯惯导解算
MINS_Ret=zeros(floor(len/2),10);kk=1;t=0;
MERR.eb=[0;0;0]*dph;
MERR.web=[0;0;0]*dpsh;
MERR.db=[0;0;0]*ug;
MERR.wdb=[0;0;0]*ugpsHz;
qnb=a2qua(att0);
vn=vn0;
posMINS=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [qnb,vn,posMINS]=my_insupdate(qnb,vn,posMINS,wm1,vm1,ts);
    MINS_Ret(kk,:)=[q2att(qnb);vn;posMINS;t];kk=kk+1;
    
    if mod(t,100)<nts,disp(fix(t));end
end
%% 子惯刚体纯惯导解算
SINS_Ret_L=zeros(floor(len/2),10);kk=1;t=0;
SERR.eb=[0;0;0]*dph;
SERR.web=[0;0;0]*dpsh;
SERR.db=[0;0;0]*ug;
SERR.wdb=[0;0;0]*ugpsHz;
qnb_L=a2qua(att0);
vn_L=vn0;
posSINS_L=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(SINS{1,1}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,1}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_L,vn_L,posSINS_L]=my_insupdate(qnb_L,vn_L,posSINS_L,wm1,vm1,ts);
    SINS_Ret_L(kk,:)=[q2att(qnb_L);vn_L;posSINS_L;t];kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
%% 子惯挠曲纯惯导解算
SINS_Ret_u=zeros(floor(len/2),10);kk=1;t=0;
qnb_u=a2qua(att0);
vn_u=vn0;
posSINS_u=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(SINS{1,2}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,2}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_u,vn_u,posSINS_u]=my_insupdate(qnb_u,vn_u,posSINS_u,wm1,vm1,ts);
    SINS_Ret_u(kk,:)=[q2att(qnb_u);vn_L;posSINS_u;t];kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
%% 绘图
delpos=deltapos(TR.pos);
delposMINS=deltapos(MINS_Ret(:,7:9));
delposSINS_L=deltapos(SINS_Ret_L(:,7:9));
delposSINS_u=deltapos(SINS_Ret_u(:,7:9));
%真实轨迹与主惯对比
figure 
plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
hold on;
plot3(delposMINS(:,1),delposMINS(:,2),delposMINS(:,3));
grid on;
xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');
err=delposMINS(end,:)-delpos(end,:);
disp '纬度误差 经度误差 高度误差'
disp(err)
%真实轨迹与子惯刚体对比
figure
plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
hold on;
plot3(delposSINS_L(:,1),delposSINS_L(:,2),delposSINS_L(:,3));
grid on;
xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');
err=delposSINS_L(end,:)-delpos(end,:);
disp '纬度误差 经度误差 高度误差'
disp(err)
%真实轨迹与子惯挠曲对比
figure
plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
hold on;
plot3(delposSINS_u(:,1),delposSINS_u(:,2),delposSINS_u(:,3));
grid on;
xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');
err=delposSINS_u(end,:)-delpos(end,:);
disp '纬度误差 经度误差 高度误差'
disp(err)
%子惯刚体与子惯挠曲对比
figure
plot3(delposSINS_L(:,1),delposSINS_L(:,2),delposSINS_L(:,3),'LineWidth',5);
hold on;
plot3(delposSINS_u(:,1),delposSINS_u(:,2),delposSINS_u(:,3));
grid on;
xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');
err=delposSINS_u(end,:)-delposSINS_L(end,:);
disp '纬度误差 经度误差 高度误差'
disp(err)
%真实轨迹与主惯误差
figure 
subplot(311)
plot(MINS_Ret(:,10),delpos(1:2:end,1)-delposMINS(:,1));ylabel('纬度距离/m')
title("真实轨迹与主惯纬度距离误差");grid on;
subplot(312)
plot(MINS_Ret(:,10),delpos(1:2:end,2)-delposMINS(:,2));ylabel('经度距离/m')
title("真实轨迹与主惯经度距离误差");grid on;
subplot(313)
plot(MINS_Ret(:,10),delpos(1:2:end,3)-delposMINS(:,3));ylabel('高/m');
grid on;
title("真实轨迹与主惯高度误差");
% %真实轨迹与子惯刚体误差
% figure 
% subplot(311)
% plot(MINS_Ret(:,10),delpos(1:2:end,1)-delposSINS_L(:,1));ylabel('纬度距离/m')
% title("真实轨迹与子惯刚体纬度距离误差");grid on;
% subplot(312)
% plot(MINS_Ret(:,10),delpos(1:2:end,2)-delposSINS_L(:,2));ylabel('经度距离/m')
% title("真实轨迹与子惯刚体经度距离误差");grid on;
% subplot(313)
% plot(MINS_Ret(:,10),delpos(1:2:end,3)-delposSINS_L(:,3));ylabel('高/m');
% grid on;
% title("真实轨迹与子惯刚体高度误差");
% %真实轨迹与子惯挠曲误差
% figure 
% subplot(311)
% plot(MINS_Ret(:,10),delpos(1:2:end,1)-delposSINS_u(:,1));ylabel('纬度距离/m')
% title("真实轨迹与子惯挠曲纬度距离误差");grid on;
% subplot(312)
% plot(MINS_Ret(:,10),delpos(1:2:end,2)-delposSINS_u(:,2));ylabel('经度距离/m')
% title("真实轨迹与子惯挠曲经度距离误差");grid on;
% subplot(313)
% plot(MINS_Ret(:,10),delpos(1:2:end,3)-delposSINS_u(:,3));ylabel('高/m');
% grid on;
% title("真实轨迹与子惯挠曲高度误差");
%主惯轨迹与子惯刚体误差
figure 
subplot(311)
plot(MINS_Ret(:,10),delposMINS(:,1)-delposSINS_L(:,1));ylabel('纬度距离/m')
title("主惯轨迹与子惯刚体纬度距离误差");grid on;
subplot(312)
plot(MINS_Ret(:,10),delposMINS(:,2)-delposSINS_L(:,2));ylabel('经度距离/m')
title("主惯轨迹与子惯刚体经度距离误差");grid on;
subplot(313)
plot(MINS_Ret(:,10),delposMINS(:,3)-delposSINS_L(:,3));ylabel('高/m');
grid on;
title("主惯轨迹与子惯刚体高度误差");
%真实轨迹与子惯挠曲误差
figure 
subplot(311)
plot(MINS_Ret(:,10),delposMINS(:,1)-delposSINS_u(:,1));ylabel('纬度距离/m')
title("主惯轨迹与子惯挠曲纬度距离误差");grid on;
subplot(312)
plot(MINS_Ret(:,10),delposMINS(:,2)-delposSINS_u(:,2));ylabel('经度距离/m')
title("主惯轨迹与子惯挠曲经度距离误差");grid on;
subplot(313)
plot(MINS_Ret(:,10),delposMINS(:,3)-delposSINS_u(:,3));ylabel('高/m');
grid on;
title("主惯轨迹与子惯挠曲高度误差");
%子惯刚体与子惯挠曲误差
figure 
subplot(311)
plot(MINS_Ret(:,10),delposSINS_L(:,1)-delposSINS_u(:,1));ylabel('纬度距离/m')
title("子惯刚体与子惯挠曲纬度距离误差");grid on;
subplot(312)
plot(MINS_Ret(:,10),delposSINS_L(:,2)-delposSINS_u(:,2));ylabel('经度距离/m')
title("子惯刚体与子惯挠曲经度距离误差");grid on;
subplot(313)
plot(MINS_Ret(:,10),delposSINS_L(:,3)-delposSINS_u(:,3));ylabel('高/m');
grid on;
title("子惯刚体与子惯挠曲高度误差");