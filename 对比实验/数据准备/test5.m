clc
clear
close all
gvar;    % 加载全局变量
ts = 0.1;
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;100];
%     俯仰角速率 横滚角速率 方位角速率 纵向加速度 持续时间
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
wat(:,1:3) = wat(:,1:3)*arcdeg/60;  % ->deg/s
[att, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
[wm, vm] = av2imu(att, vn, pos, ts);
tt = (0:length(att)-1)'*ts;
% 轨迹作图
msplot(221, tt, att/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
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
figure(3)
delpos=deltapos(pos);
plot3(delpos(:,1),delpos(:,2),delpos(:,3));
grid on;
xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');

qnb=a2qua(att0);
vn=vn0;
posSINS=pos0;

% posSINS(3)=110;

len=length(vm)-1;
nn=2;nts=nn*ts;
SINS=zeros(floor(len/2),10);kk=1;t=0;
% eb=[0.01;0.01;0.01]*dph;
% web=[0.001;0.001;0.001]*dpsh;
% db=[10;10;10]*ug;
% wdb=[0;0;0]*ugpsHz;
eb=[0;0;0]*dph;
web=[0;0;0]*dpsh;
db=[0;0;0]*ug;
wdb=[0;0;0]*ugpsHz;
sensor1=zeros(2,6,len/2);
sensor=zeros(2,6,len/2);
num=1;
for k=1:nn:len
    t=t+nts;
    [wm1,vm1]=imuadderr(wm(k:k+1,:),vm(k:k+1,:), eb, web, db, wdb, ts);
    sensor1(:,:,num)=[wm1,vm1];sensor(:,:,num)=[wm1,vm1];
    num=num+1;
    [qnb,vn,posSINS]=insupdate(qnb,vn,posSINS,wm1,vm1,ts);
    SINS(kk,:)=[q2att(qnb);vn;posSINS;t];kk=kk+1;
    
    if mod(t,100)<nts,disp(fix(t));end
end
%%
figure(4)
delpos=deltapos(pos);
plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
hold on;
delposSINS=deltapos(SINS(:,7:9));
plot3(delposSINS(:,1),delposSINS(:,2),delposSINS(:,3));
grid on;
xlabel('纬度距离/m'),ylabel('经度距离/m'),zlabel('高/m');
err=delposSINS(end,:)-delpos(end,:);
disp '纬度误差 经度误差 高度误差'
disp(err)