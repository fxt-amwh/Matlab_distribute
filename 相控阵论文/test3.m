%% 相控阵位置误差加入后峰值衰减
% 版本时间：2020/12/19 
% 设计初衷：机翼上弯曲变形对蒙皮天线的影响
clc
clear
close all
global lambda
global M
global N
sita=linspace(-pi/2,pi/2,2000);
phi=-pi/2;
cm=0.01;
M=33;%x方向阵元数量
N=9;%y方向阵于M
lambda=3*cm;%波长
dx=lambda*2;%x方向阵元距离
dy=lambda*2;%y方向阵元距离
xm=floor(M/2)*dx;%x方向最大距离
ym=floor(N/2)*dy;%y方向最大距离
f=my_getDirPt(sita,phi,0,0,0);%无位置误差
flag=1;%控制仿真模式 1：x方向的弯曲变形 2：弯曲变形+随机二维位姿误差 
SigmaMax=0.1*lambda;%控制仿真位置误差最大值
SitaMax=1;%控制仿真俯仰误差最大值
PhiMax=1;%控制仿真方向误差最大值
simLen=100;%蒙特卡洛仿真次数
Zmax=1*(0:0.025:0.5)*lambda;%控制变形指数
% Zmax=1*(0:0.02:0.4)*lambda;%控制变形指数
sigmaList=Zmax/max(Zmax)*SigmaMax;
sigmasitaList=Zmax/max(Zmax)*SitaMax*pi/180;
sigmaphiList=Zmax/max(Zmax)*PhiMax*pi/180;
if flag<2
    dGList=zeros(1,length(Zmax));
    feList=cell(size(Zmax));
    for loopi=1:length(Zmax)
        feList{loopi}=zeros(size(f));
    end
else
    dGList=zeros(length(Zmax),length(sigmasitaList));
    feList=cell(length(Zmax),length(sigmasitaList));
    for loopi=1:length(Zmax)
        for loopj=1:length(sigmasitaList)
            feList{loopi,loopj}=zeros(size(f));
        end
    end
end

tic
for loopi=1:length(Zmax)
    if simLen>1
         fprintf('仿真中');
    end
    for simi=1:simLen
        switch flag
            case 1
                ddx=0*sigmaList(loopi)*randn(N,M);
                ddy=0*sigmaList(loopi)*randn(N,M);
                ddz=0*sigmaList(loopi)*randn(N,M);
                xdir=(-floor(M/2):1:floor(M/2))*dx;
                dz=my_fixdef12(xdir,xm,ym,Zmax(loopi),1);
                nowfe=my_getDirPt(sita,phi,ddx,ddy,dz+ddz);
                feList{loopi}=feList{loopi}+nowfe./simLen;%位置误差
                dGList(loopi)=dGList(loopi)-20*log10(max(nowfe)/max(f));
            case 2
                for loopj=1:length(Zmax)
                    ddx=sigmaList(loopi)*randn(N,M);
                    ddy=sigmaList(loopi)*randn(N,M);
                    ddz=sigmaList(loopi)*randn(N,M);
                    xdir=(-floor(M/2):1:floor(M/2))*dx;
                    dz=my_fixdef(xdir,xm,ym,Zmax(loopj),1,ddx,ddy);
                    nowfe=my_getDirPt(sita,phi,ddx,ddy,dz+ddz);
                    feList{loopi,loopj}=feList{loopi,loopj}+nowfe./simLen;%位置、姿态误差
                    dGList(loopi,loopj)=dGList(loopi,loopj)-20*log10(max(nowfe)/max(f));
                end
        end
        if simLen>1
            if mod(simi,ceil(simLen/10))==0
                fprintf('.');
            end
        end
    end
    if flag<2
        dGList(loopi)=dGList(loopi)/simLen;
    else
        for loopj=1:length(Zmax)
            dGList(loopi,loopj)=dGList(loopi,loopj)/simLen;
        end
    end
    fprintf(' ');
    switch flag
        case 1 
            fprintf('变形指数%f时峰值衰减：%f dB\n'...
            ,Zmax(loopi)/lambda,dGList(loopi))
        case 2
            fprintf('位置误差标准差：%f，变形指数：%f时峰值衰减：%f dB\n'...
            ,sigmaList(loopi)/lambda,Zmax(loopi)/lambda,dGList(loopi,loopi))
    end
end
toc
%%
figure
if flag<2
    plot(sita,f(1,:)/max(f),'y','LineWidth',3);
    hold on;
    plot(sita,feList{end,end}/max(f),'r');
else
    plot(sita,f(1,:)/max(f),'y','LineWidth',3);
    hold on;
    plot(sita,feList{end,end}/max(f),'r');
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
        %拟合
        [xi,yi,fuc]=my_fit(Zmax/lambda,dGList);
        %绘图
        figure
        plot(Zmax/lambda,dGList,'*');
        xlabel('变形指数/\lambda','Fontsize',15);
        hold on;plot(xi,yi);
        ylabel('增益损失/dB','Fontsize',15);
        title('弯曲变形增益损失图','Fontsize',15);
        legend('仿真','拟合');
    case 2
        figure
        x=repmat(sigmaList/lambda,length(Zmax),1);
        y=repmat(Zmax'/lambda,1,length(sigmaList));
        z=dGList';
        surf(x,y,z);
        xlabel('位置随机误差均方差/\lambda','Fontsize',15);
        ylabel('变形指数/\lambda','Fontsize',15);
        zlabel('增益损失/dB','Fontsize',15);
        title('弯曲变形+随机误差增益损失图','Fontsize',15);
end
if flag<2
    disp('拟合函数:')
    fuc
end
%%
figure
ydir=(-floor(N/2):1:floor(N/2))*dy;
[MList,NList]=meshgrid(xdir,ydir);
switch mod(flag,3)
    case 1
        surf((MList+ddx)/cm,(NList+ddy)/cm,(dz+ddz)/cm);
        title('弯曲变形','Fontsize',15);
end
xlabel('x/cm','Fontsize',15);ylabel('y/cm','Fontsize',15);zlabel('z/cm','Fontsize',15);
