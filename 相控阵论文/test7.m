%% 相控阵位置误差加入后峰值衰减 14所建议后修改
% 版本时间：2020/01/21 
% 设计初衷：机翼上弯曲变形对蒙皮天线的影响 封装 2m处中心安装蒙皮 精细地解算位移
% 论文图，flag 1 200 s - 500 s 的补偿结果 - 最大值 -均值
% 论文图，flag 2 200 s - 500 s 的补偿结果 - 每一个采样时刻对一个随机安装误差角进行多次蒙特卡洛仿真先求蒙特卡洛的均值
% 再在 200 s - 500 s 求最大值 均值
clc
clear
close all
sita=linspace(-pi/2,pi/2,2000);
phi=0;
cm=0.01;
M=33;%x方向阵元数量
N=21;%y方向阵于M
lambda=5*cm;%波长
dx=lambda*0.6;%x方向阵元距离
dy=lambda*0.6;%y方向阵元距离
xm=floor(M/2)*dx;%x方向最大距离
ym=floor(N/2)*dy;%y方向最大距离
RO.base.dx=lambda*0.6;%x方向阵元距离
RO.base.dy=lambda*0.6;%y方向阵元距离
RO.base.phi=phi;
RO.base.sita=sita;
RO.base.lambda=lambda;%波长
RO.base.M=M;%x方向阵元个数
RO.base.N=N;%y方向阵元个数
RO.base.sita0=0;%波束控制指向俯仰角
RO.base.phi0=0;%波束控制指向方位角
RO.base.ddsita=0;%俯仰角误差矩阵
RO.base.ddphi=0;%方位角误差矩阵
RO.base.ddx=zeros(N,M);%x方向位置误差矩阵
RO.base.ddy=zeros(N,M);%y方向位置误差矩阵
RO.base.ddz=zeros(N,M);%z方向位置误差矩阵

RO.f=my_getDirPtFoc(RO.base);%无位置误差
flag=1;%控制仿真模式 1：x方向的弯曲变形 2：弯曲变形+随机二维位姿误差 
SigmaMax=0.1*lambda;%控制仿真位置误差最大值
SitaMax=1;%控制仿真俯仰误差最大值
PhiMax=1;%控制仿真方向误差最大值
simLen=1;%蒙特卡洛仿真次数
Zmax=1*(0:0.025:0.5)*lambda;%控制变形指数
% umax=Zmax/max(Zmax)*5*pi/180;%控制变形指数
ts=0.1;%采样时间
nts=2*ts;
samp.begin=200/nts;
samp.end=500/nts;
samp.dir=10*2;
load NSINSFR3;
% umax=NSINSFR.attnewall_true(2,1:2:7000)/2;%单位机翼处挠曲角真实值 载荷模型
% umax_measure=NSINSFR.attnewall(2,1:2:7000)/2;%单位机翼处挠曲角量测值 载荷模型
% umax=NSINSFR.attnewall_true(2,:)/2;%单位机翼处挠曲角真实值 载荷模型
% umax_measure=NSINSFR.attnewall(2,:)/2;%单位机翼处挠曲角量测值 载荷模型
umax=NSINSFR.attnewall_true(2,samp.begin:samp.dir:samp.end)/2;%单位机翼处挠曲角真实值 载荷模型
umax_measure=NSINSFR.attnewall(2,samp.begin:samp.dir:samp.end)/2;%单位机翼处挠曲角量测值 载荷模型
sigmaList=Zmax/max(Zmax)*SigmaMax;
% sigmasitaList=Zmax/max(Zmax)*SitaMax*pi/180;
% sigmaphiList=Zmax/max(Zmax)*PhiMax*pi/180;

if flag==1
    len_loopi=length(umax);
elseif flag==2
    len_loopi=length(sigmaList);
    len_loopj=length(umax);
end
if flag<2
    simLen=1;
    dGList=zeros(1,len_loopi);
    RO.feList=cell(len_loopi,1);
    for loopi=1:len_loopi
        RO.feList{loopi}=zeros(size(RO.f));
    end
else
    simResult=zeros(len_loopi,2);
    dGList=zeros(len_loopi,len_loopj);
    RO.feList=cell(len_loopi,len_loopj);
    for loopi=1:len_loopi
        for loopj=1:len_loopj
            RO.feList{loopi,loopj}=zeros(size(RO.f));
        end
    end
end

tic

for loopi=1:len_loopi
    if flag==1

    elseif flag==2
        fprintf('仿真...%6.1f %%%6.1f %%',0,0);
    end
    if simLen>1
         fprintf('%6.1f %%',0);
    end
    for simi=1:simLen
        switch flag
            case 1
                ddx=0*0.001*randn(N,M);
                ddy=0*0.001*randn(N,M);
                ddz=0*0.001*randn(N,M);
                xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%蒙皮中心安装在 2m 处
                [d_x,dz]=my_fixdefload(xdir+ddx,1,umax(loopi),1);%机翼载荷挠曲模型标准量 精准模型1
                [dx_measure,dz_measure]=my_fixdefload(xdir+ddx,1,umax_measure(loopi),1);%机翼载荷挠曲模型补偿量
                RO.base.ddx=d_x;%x方向位置误差矩阵
                RO.base.ddy=ddy;%y方向位置误差矩阵
                RO.base.ddz=dz+ddz;%z方向位置误差矩阵
                RO.base.ddx=d_x-dx_measure;%x方向位置误差矩阵
                RO.base.ddz=dz+ddz-dz_measure;%z方向位置误差矩阵 补偿
                
                nowfe=my_getDirPtFoc(RO.base);
                RO.feList{loopi}=RO.feList{loopi}+nowfe./simLen;%位置误差
                dGList(loopi)=dGList(loopi)-20*log10(max(nowfe)/max(RO.f));
            case 2
                for loopj=1:len_loopj
                    ddx=sigmaList(loopi)*randn(N,M);
                    ddy=sigmaList(loopi)*randn(N,M);
                    ddz=sigmaList(loopi)*randn(N,M);
                    xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%蒙皮中心安装在 2m 处
                    [d_x,dz]=my_fixdefload(xdir+ddx,1,umax(loopj),1);%机翼载荷挠曲模型标准量 精准模型1
                    [dx_measure,dz_measure]=my_fixdefload(xdir+ddx,1,umax_measure(loopj),1);%机翼载荷挠曲模型补偿量
                    RO.base.ddx=d_x;%x方向位置误差矩阵
                    RO.base.ddy=ddy;%y方向位置误差矩阵
                    RO.base.ddz=dz+ddz;%z方向位置误差矩阵
                    RO.base.ddx=d_x-dx_measure;%x方向位置误差矩阵
                    RO.base.ddz=dz+ddz-dz_measure;%z方向位置误差矩阵 补偿
                
%                 RO.base.ddz=dz+ddz-dz_measure;%z方向位置误差矩阵 补偿
                
                    nowfe=my_getDirPtFoc(RO.base);
                    RO.feList{loopi,loopj}=RO.feList{loopi,loopj}+nowfe./simLen;%位置、姿态误差
                    dGList(loopi,loopj)=dGList(loopi,loopj)-20*log10(max(nowfe)/max(RO.f));
                    if simLen==1
                        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6.1f %%%6.1f %%', ((loopi-1)*len_loopj+loopj)/(len_loopi*len_loopj)*100, loopj/len_loopj*100);	
                    else
                        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6.1f %%%6.1f %%%6.1f %%', ((loopi-1)*simLen*len_loopj+(simi-1)*len_loopj+loopj)/(len_loopi*simLen*len_loopj)*100,simi/simLen*100,loopj/len_loopj*100);	
                    end
                end
        end
        if simLen>1&&flag==1
            if mod(simi,ceil(simLen/10))==0
                fprintf('.');
            end
        end
    end
    
    if flag<2
        dGList(loopi)=dGList(loopi)/simLen;
        
    else
        for loopj=1:len_loopj
            dGList(loopi,loopj)=dGList(loopi,loopj)/simLen;
            simResult(loopi,:)=[max(dGList(loopi,:)),mean(dGList(loopi,:))];
        end
        fprintf("\n");
    end
    switch flag
        case 1 
            fprintf('%d挠曲指数%f时峰值衰减：%f dB\n'...
            ,loopi,5*umax(loopi)*180/pi,dGList(loopi))
        case 2
            fprintf('位置误差标准差：%f，峰值衰减最大值：%f dB 均值：%f dB\n'...
            ,sigmaList(loopi)/lambda,simResult(loopi,1),simResult(loopi,2))
    end
end
toc
%%
figure
if flag<2
    plot(sita,RO.f(1,:)/max(RO.f),'y','LineWidth',3);
    hold on;
    plot(sita,RO.feList{end,end}/max(RO.f),'r');
else
    plot(sita,RO.f(1,:)/max(RO.f),'y','LineWidth',3);
    hold on;
    plot(sita,RO.feList{end,end}/max(RO.f),'r');
end
grid on;
legend('无误差','随机误差');
xlabel('theta/radian');
ylabel('amplitude');
title('归一化方向图phi=0')
%% 三维方向图

%%
switch flag
    case 1
        %绘图
        figure
        plot(linspace(samp.begin*nts,samp.end*nts,len_loopi),dGList,'LineWidth',2);
        xlabel('时间/s','Fontsize',15);
        ylabel('增益损失/dB','Fontsize',15);
%         title('机翼挠曲变形未补偿增益损失图','Fontsize',15);
%         title('机翼挠曲变形补偿后增益损失图(\sigma=0.001mm)','Fontsize',15);
        title('机翼挠曲变形补偿后增益损失图','Fontsize',15);
        fprintf("最大增益损失%fdB，平均增益损失%fdB\n",max(dGList),mean(dGList));
        %绘图
        figure
        plot(linspace(samp.begin*nts,samp.end*nts,len_loopi),dGList,'LineWidth',2);
        hold on;
        plot(linspace(samp.begin*nts,samp.end*nts,len_loopi),max(dGList)*ones(1,len_loopi),"-.",'LineWidth',2);
        hold on;
        plot(linspace(samp.begin*nts,samp.end*nts,len_loopi),mean(dGList)*ones(1,len_loopi),"--",'LineWidth',2);
        xlabel('时间/s','Fontsize',15);
        ylabel('增益损失/dB','Fontsize',15);
        title('机翼挠曲变形未补偿增益损失图','Fontsize',15);
%         title('机翼挠曲变形补偿后增益损失图(\sigma=0.001mm)','Fontsize',15);
        legend("补偿后","最大值","均值");
    case 2
        figure
        x=repmat(sigmaList/lambda,len_loopj,1);
        y=repmat(linspace(samp.begin*nts,samp.end*nts,len_loopj)',1,len_loopi);
        z=dGList';
        surf(x,y,z);
        xlabel('位置随机误差均方差/\lambda','Fontsize',15);
        ylabel('时间/s','Fontsize',15);
        zlabel('增益损失/dB','Fontsize',15);
        title('机翼挠曲变形+随机误差增益损失图','Fontsize',15); 
        figure
        plot(sigmaList/lambda,simResult(:,1),'LineWidth',2);
        hold on;
        plot(sigmaList/lambda,simResult(:,2),'LineWidth',2);
        xlabel('安装位置误差均方差/\lambda','Fontsize',15);
        ylabel('增益损失/dB','Fontsize',15);
        title('相位补效果','Fontsize',15);
        legend("最大值","平均值");
end
%%
figure
ydir=(-floor(N/2):1:floor(N/2))*dy;
[MList,NList]=meshgrid(xdir,ydir);
[d_x,dz]=my_fixdefload(xdir+ddx,5,1*pi/180,1);%机翼载荷挠曲模型标准量 精准模型1
switch mod(flag,3)
    case 1
        subplot(211)
        surf((MList+d_x+0*ddx)/cm,(NList+0*ddy)/cm,(dz+0*ddz)/cm);
%         title('机翼挠曲变形图','Fontsize',15);
        title('挠曲指数1°时天线阵元位置','Fontsize',15);
        ylim([-50,50]);zlim([0,10]);
        colorbar;
        subplot(212)
        surf((MList+d_x+ddx)/cm,(NList+0*ddy)/cm,(dz+ddz)/cm);
%         title('机翼挠曲变形图','Fontsize',15);
        title('挠曲指数1°时天线阵元位置(\sigma=0.001mm)','Fontsize',15);
        ylim([-50,50]);zlim([0,10]);
        colorbar;
end
xlabel('x/cm','Fontsize',15);ylabel('y/cm','Fontsize',15);zlabel('z/cm','Fontsize',15);
