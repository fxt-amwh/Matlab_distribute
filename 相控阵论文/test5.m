%% 相控阵位置误差加入后峰值衰减
% 版本时间：2020/12/20 
% 设计初衷：机翼上弯曲变形对蒙皮天线的影响 封装 2m处中心安装蒙皮 精细地解算位移
clc
clear
close all
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
simLen=100;%蒙特卡洛仿真次数
Zmax=1*(0:0.025:0.5)*lambda;%控制变形指数
umax=Zmax/max(Zmax)*0.5*pi/180;%控制变形指数
sigmaList=Zmax/max(Zmax)*SigmaMax;
sigmasitaList=Zmax/max(Zmax)*SitaMax*pi/180;
sigmaphiList=Zmax/max(Zmax)*PhiMax*pi/180;
len=length(umax);
if flag<2
    simLen=1;
    dGList=zeros(1,len);
    RO.feList=cell(len,1);
    for loopi=1:len
        RO.feList{loopi}=zeros(size(RO.f));
    end
else
    dGList=zeros(len,len);
    RO.feList=cell(len,len);
    for loopi=1:len
        for loopj=1:len
            RO.feList{loopi,loopj}=zeros(size(RO.f));
        end
    end
end

tic
for loopi=1:len
    if simLen>1
         fprintf('仿真中');
    end
    for simi=1:simLen
        switch flag
            case 1
                ddx=0*sigmaList(loopi)*randn(N,M);
                ddy=0*sigmaList(loopi)*randn(N,M);
                ddz=0*sigmaList(loopi)*randn(N,M);
                xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%蒙皮中心安装在 2m 处
                [d_x,dz]=my_fixdefload(xdir,5,umax(loopi),1);%机翼载荷挠曲模型标准量 精准模型1
                [dx_measure,dz_measure]=my_fixdefload(xdir,5,umax(loopi)+0.1*pi/180,1);%机翼载荷挠曲模型补偿量
                RO.base.ddx=d_x+ddx-dx_measure;%x方向位置误差矩阵
                RO.base.ddy=ddy;%y方向位置误差矩阵
                RO.base.ddz=dz+ddz;%z方向位置误差矩阵
                
                RO.base.ddz=dz+ddz-dz_measure;%z方向位置误差矩阵 补偿
                
                nowfe=my_getDirPtFoc(RO.base);
                RO.feList{loopi}=RO.feList{loopi}+nowfe./simLen;%位置误差
                dGList(loopi)=dGList(loopi)-20*log10(max(nowfe)/max(RO.f));
            case 2
                for loopj=1:len
                    ddx=sigmaList(loopi)*randn(N,M);
                    ddy=sigmaList(loopi)*randn(N,M);
                    ddz=sigmaList(loopi)*randn(N,M);
                    xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%蒙皮中心安装在 2m 处
                    [d_x,dz]=my_fixdefload(xdir,5,umax(loopj),1);%机翼载荷挠曲模型标准量 精准模型1
                    dz_measure=my_fixdefload(xdir,5,umax(loopj)-0.1*pi/180);%机翼载荷挠曲模型补偿量
                    RO.base.ddx=d_x+ddx;%x方向位置误差矩阵
                    RO.base.ddy=ddy;%y方向位置误差矩阵
                    RO.base.ddz=dz+ddz;%z方向位置误差矩阵
                
%                 RO.base.ddz=dz+ddz-dz_measure;%z方向位置误差矩阵 补偿
                
                    nowfe=my_getDirPtFoc(RO.base);
                    RO.feList{loopi,loopj}=RO.feList{loopi,loopj}+nowfe./simLen;%位置、姿态误差
                    dGList(loopi,loopj)=dGList(loopi,loopj)-20*log10(max(nowfe)/max(RO.f));
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
        for loopj=1:len
            dGList(loopi,loopj)=dGList(loopi,loopj)/simLen;
        end
    end
    fprintf(' ');
    switch flag
        case 1 
            fprintf('挠曲角%f时峰值衰减：%f dB\n'...
            ,umax(loopi)*180/pi,dGList(loopi))
        case 2
            fprintf('位置误差标准差：%f，挠曲角：%f时峰值衰减：%f dB\n'...
            ,sigmaList(loopi)/lambda,umax(loopi)*180/pi,dGList(loopi,loopi))
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
        %拟合
        [xi,yi,fuc]=my_fit(umax*180/pi,dGList);
        %绘图
        figure
        plot(umax*180/pi,dGList,'*');
        xlabel('挠曲指数/°','Fontsize',15);
        hold on;plot(xi,yi);
        ylabel('增益损失/dB','Fontsize',15);
        title('机翼挠曲变形0.1°补偿后增益损失图','Fontsize',15);
        legend('仿真','拟合');
    case 2
        figure
        x=repmat(sigmaList/lambda,len,1);
        y=repmat(umax'*180/pi,1,len);
        z=dGList';
        surf(x,y,z);
        xlabel('位置随机误差均方差/\lambda','Fontsize',15);
        ylabel('挠曲指数/°','Fontsize',15);
        zlabel('增益损失/dB','Fontsize',15);
        title('机翼挠曲变形+随机误差增益损失图','Fontsize',15);
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
%%
figure
ddx=0.001*randn(N,M)*0;
ddy=0.001*randn(N,M)*0;
ddz=0.001*randn(N,M)*0;
ydir=(-floor(N/2):1:floor(N/2))*dy;
[MList,NList]=meshgrid(xdir,ydir);
[d_x,dz]=my_fixdefload(xdir+ddx,5,1*pi/180,1);%机翼载荷挠曲模型标准量 精准模型1
switch mod(flag,3)
    case 1
%         subplot(211) 
        plot3((MList+d_x+ddx)/cm,(NList+ddy)/cm,(dz+ddz)/cm,"r+");
        hold on;
        surf((MList+d_x+ddx)/cm,(NList+ddy)/cm,(dz+ddz)/cm);
%         title('机翼挠曲变形图','Fontsize',15);
%         title('挠曲指数1°时天线阵元位置','Fontsize',15);
        ylim([-50,50]);
        zlim([0,5]);
        colorbar;
        legend("阵元")
        grid on;
%         subplot(212)
%         surf((MList+d_x+ddx)/cm,(NList+0*ddy)/cm,(dz+ddz)/cm);
% %         title('机翼挠曲变形图','Fontsize',15);
%         title('挠曲指数1°时天线阵元位置(\sigma=0.001mm)','Fontsize',15);
%         ylim([-50,50]);zlim([0,10]);
%         colorbar;
end
xlabel('x/cm','Fontsize',15);ylabel('y/cm','Fontsize',15);zlabel('z/cm','Fontsize',15);
