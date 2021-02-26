%% 一主+多子（高度模块化）for汇报 蒙特卡洛
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
simLen=1;
simListLen=4;
gvar;    % 加载全局变量
ts = 0.1;%采样时间
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
%         ];    %静止
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

for i=1:simListLen
    nameList(i)=sprintf("SINS%d",i);
end
% 子惯运动数据设计
fprintf('子惯仿真...'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;
Sinf.Rlist=[0.5*ones(1,simListLen);
            0*ones(1,simListLen);
            0*ones(1,simListLen)];
Sinf.ulist=2.5*ones(1,simListLen)*arcdeg;%子惯处挠曲角
Sinf.flist=1*ones(1,simListLen)*0.01;%挠曲频率
Sinf.aerrlist=repmat([0;0;0]*arcdeg,1,simListLen);%子惯安装误差角
Sinf.wim0=MINS.wim0;
SinfCell{1,1}=Sinf;
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
    atterrList(:,:,m)=repmat([0.01;0.01;0.01]*arcdeg,1,simListLen);%子惯安装误差角
end
%% 子惯反解算数据
SINS=my_invRIpackN(MINS,Smove);
fprintf('子惯仿真完成！\n'); 
%% 子惯误差转换为元组
SERR=my_getSERRN(SERRList);
%% 蒙特卡洛
result=zeros(18,N*simLen,M);
resultF=zeros(18,simLen);
resultR=zeros(3,3599);
resultAtt=zeros(3,3599);
tic
for loopi=1:simLen
    %% 相对导航
    disp(loopi);
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
            [SINSFR,Filter,NSINSFR]=my_getFResultLoopN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,flag);%并行集中版
        else
            [SINSFR,Filter,~]=my_getFResultLoopN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,flag);%并行集中版
        end
    else
        flag.select=2;
        % [SINSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,2);%单个子惯分离版
        if flag.EnFusion
            [SINSFR,Filter,NSINSFR]=my_getFResultLoopN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,flag);%并行集中版
        else
            [SINSFR,Filter,~]=my_getFResultLoopN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,flag);%并行集中版
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
%             disp(nameList(m,n));
            errlen=1/(Smove{m,n}.f*nts);
            Rerr=my_getSErr(SINSFR{m,n}.Rall(1:end,:)'-SINS{m,n}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
            Atterr=my_getSErr(SINSFR{m,n}.attall(:,1:end)-SINS{m,n}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
%             fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%                 Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
            result(7:9,n+N*(loopi-1),m)=Rerr;result(10:12,n+N*(loopi-1),m)=Atterr;
            errlen=1;
            beg=ceil(200/nts);
            End=ceil(500/nts);
            Rerr=my_getSErr(SINSFR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
            Rserr=my_getVErr(SINSFR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),beg,End,2)*1000;        %位置误差单位mm
            Atterr=my_getSErr(SINSFR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
            Attserr=my_getVErr(SINSFR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),beg,End,2)/arcdeg*60;
%             fprintf('滤波终值,位置误差%f %f %f（mm）,位置误差方差%f %f %f(mm)，姿态误差%f %f %f（分）,姿态均方差%f %f %f(分)\n',...
%                 Rerr(1),Rerr(2),Rerr(3),Rserr(1),Rserr(2),Rserr(3),Atterr(1),Atterr(2),Atterr(3),Attserr(1),Attserr(2),Attserr(3));
            result(1:3,n+N*(loopi-1),m)=Rerr;result(4:6,n+N*(loopi-1),m)=Atterr;result(13:15,n+N*(loopi-1),m)=Rserr;result(16:18,n+N*(loopi-1),m)=Attserr;
            if flag.EnReli
                errlen=1/(Smove{m,n}.f*nts);
                Rerr=my_getSErr(SINSR{m,n}.Rall(1:end,:)'-SINS{m,n}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
                Atterr=my_getSErr(SINSR{m,n}.attall(:,1:end)-SINS{m,n}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
%                 fprintf('纯惯统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%                     Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
                errlen=1;
                Rerr=my_getSErr(SINSR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
                Atterr=my_getSErr(SINSR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
%                 fprintf('纯惯终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%                     Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
            end
        end
    end

    if flag.EnFusion
%         disp('融合结果');
    %     errlen=1/(Smove{1,1}.f*nts);
    %     Rerr=my_getSErr(NSINSFR.Rall(1:end,:)'-SINS{1,1}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
    %     Atterr=my_getSErr(NSINSFR.attall(:,1:end)-SINS{1,1}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
    %     fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    %         Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
    %     errlen=1;
    %     Rerr=my_getSErr(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
    %     Atterr=my_getSErr(NSINSFR.attall(:,1:end-1)-SINS{1,1}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
    %     fprintf('滤波终值,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
    %         Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
        errlen=1/(Smove{1,1}.f*nts);
        Rerr=my_getSErr(NSINSFR.Rall(1:end,:)'-SINS{1,1}.R(:,1:2:end),errlen,2)*1000;        %位置误差单位mm
        Atterr=my_getSErr(NSINSFR.attall(:,1:end)-SINS{1,1}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %姿态误差单位分
%         fprintf('滤波统计,位置误差%f %f %f（mm），姿态误差%f %f %f（分）\n',...
%             Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
        resultF(7:9,loopi)=Rerr;resultF(10:12,loopi)=Atterr;
        errlen=1;
        beg=ceil(200/nts);
        End=ceil(500/nts);
        Rerr=my_getSErr(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),errlen,2)*1000;        %位置误差单位mm
        Rserr=my_getVErr(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),beg,End,2)*1000;        %位置误差单位mm
        Atterr=my_getSErr(NSINSFR.attall(:,1:end-1)-SINS{1,1}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %姿态误差单位分
        Attserr=my_getVErr(NSINSFR.attall(:,1:end-1)-SINS{1,1}.atttrue(:,1:2:end-2),beg,End,2)/arcdeg*60;
%         fprintf('滤波终值,位置误差%f %f %f（mm）,位置误差方差%f %f %f(mm)，姿态误差%f %f %f（分）,姿态均方差%f %f %f(分)\n',...
%             Rerr(1),Rerr(2),Rerr(3),Rserr(1),Rserr(2),Rserr(3),Atterr(1),Atterr(2),Atterr(3),Attserr(1),Attserr(2),Attserr(3));   
        resultF(1:3,loopi)=Rerr;resultF(4:6,loopi)=Atterr;resultF(13:15,loopi)=Rserr;resultF(16:18,loopi)=Attserr;
        
        resultR=resultR+NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2);
        resultAtt=resultAtt+NSINSFR.attall(:,1:end-1)-SINS{1,1}.atttrue(:,1:2:end-2);
        sam.begin=200/nts;
        sam.end=500/nts;
        RMSE=1000*[sqrt(diag(resultR(:,sam.begin:sam.end)*resultR(:,sam.begin:sam.end)')/length(resultR(:,sam.begin:sam.end))/loopi^2);sqrt(diag(resultAtt(:,sam.begin:sam.end)*resultAtt(:,sam.begin:sam.end)')/length(resultAtt(:,sam.begin:sam.end))/loopi^2)];
        fprintf("RMSE:(mm,mrad)");
        disp(RMSE')
    end
end
%%
sam.begin=200/nts;
sam.end=500/nts;
RMSE=1000*[sqrt(diag(resultR(:,sam.begin:sam.end)*resultR(:,sam.begin:sam.end)')/length(resultR(:,sam.begin:sam.end))/simLen^2);sqrt(diag(resultAtt(:,sam.begin:sam.end)*resultAtt(:,sam.begin:sam.end)')/length(resultAtt(:,sam.begin:sam.end))/simLen^2)];
fprintf("RMSE:(mm,mrad)");
disp(RMSE')
toc