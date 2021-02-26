%% 相控阵位置误差加入后峰值衰减
% 版本时间：2021/1/23 
% 设计初衷：机翼上弯曲变形对蒙皮天线的影响 封装 2m处中心安装蒙皮 精细地解算位移
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
RO.base.dx=lambda*0.6;%x方向阵元距离
RO.base.dy=lambda*0.6;%y方向阵元距离
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
% umax=Zmax/max(Zmax)*5*pi/180;%控制变形指数
load NSINSFR8;
% umax=NSINSFR.attnewall_true(2,1:2:7000)/2;%单位机翼处挠曲角真实值 载荷模型
% umax_measure=NSINSFR.attnewall(2,1:2:7000)/2;%单位机翼处挠曲角量测值 载荷模型
dip=1;
umax=NSINSFR.attnewall_true(2,1:dip:end)/2;%单位机翼处挠曲角真实值 载荷模型
umax_measure=NSINSFR.attnewall(2,1:dip:end)/2;%单位机翼处挠曲角量测值 载荷模型

sigmaList=umax/max(umax)*SigmaMax;
sigmasitaList=umax/max(umax)*SitaMax*pi/180;
sigmaphiList=umax/max(umax)*PhiMax*pi/180;
len=length(umax);
if flag<2
    simLen=1;
    dGList=zeros(1,len);
    dGList_err=zeros(1,len);
    sitaList=zeros(1,len);
    sitaList_err=zeros(1,len);
    RO.feList=cell(len,1);
    for loopi=1:len
        RO.feList{loopi}=zeros(size(RO.f));
    end
else
    dGList=zeros(len,len);
    dGList_err=zeros(len,len);
    sitaList=zeros(len,len);
    sitaList_err=zeros(len,len);
    RO.feList=cell(len,len);
    for loopi=1:len
        for loopj=1:len
            RO.feList{loopi,loopj}=zeros(size(RO.f));
        end
    end
end
RO_err=RO;
tic
for loopi=1:len
    if simLen>1
         fprintf('仿真中');
    end
    for simi=1:simLen
        switch flag
            case 1
                ddx=0.001*randn(N,M);
                ddy=0.001*randn(N,M);
                ddz=0.001*randn(N,M);
                xdir=0.5 + (-floor(M/2):1:floor(M/2))*dx;%蒙皮中心安装在 2m 处
                [d_x,dz]=my_fixdefload(xdir,2,2*umax(loopi),1);%机翼载荷挠曲模型标准量 精准模型1
                [dx_measure,dz_measure]=my_fixdefload(xdir,2,2*umax_measure(loopi),1);%机翼载荷挠曲模型补偿量
                
                RO.base.ddx=d_x+ddx-dx_measure;%x方向位置误差矩阵
                RO.base.ddy=ddy;%y方向位置误差矩阵
                RO.base.ddz=dz+ddz-dz_measure;%z方向位置误差矩阵 补偿
                
                RO_err.base.ddx=d_x+ddx;%x方向位置误差矩阵
                RO_err.base.ddy=ddy;%y方向位置误差矩阵
                RO_err.base.ddz=dz+ddz;%z方向位置误差矩阵
                
                nowfe=my_getDirPtFoc(RO.base);%补偿后
                nowfe_err=my_getDirPtFoc(RO_err.base);%不补偿
                RO.feList{loopi}=RO.feList{loopi}+nowfe./simLen;%位置误差
                RO_err.feList{loopi}=RO_err.feList{loopi}+nowfe_err./simLen;%位置误差
                [maxfe,maxfedir]=max(nowfe);
                [maxfe_err,maxfedir_err]=max(nowfe_err);
                dGList(loopi)=dGList(loopi)-20*log10(maxfe/max(RO.f));
                dGList_err(loopi)=dGList_err(loopi)-20*log10(maxfe_err/max(RO_err.f));
                
                sitaList(loopi)=sita(maxfedir)*180/pi;%记录当前sita误差
                sitaList_err(loopi)=sita(maxfedir_err)*180/pi;
            case 2
                for loopj=1:len
                    ddx=sigmaList(loopi)*randn(N,M);
                    ddy=sigmaList(loopi)*randn(N,M);
                    ddz=sigmaList(loopi)*randn(N,M);
                    xdir=0.5 + (-floor(M/2):1:floor(M/2))*dx;%蒙皮中心安装在 2m 处
                    [d_x,dz]=my_fixdefload(xdir,2,2*umax(loopj),1);%机翼载荷挠曲模型标准量 精准模型1
                    [dx_measure,dz_measure]=my_fixdefload(xdir,2,2*umax_measure(loopj),1);%机翼载荷挠曲模型补偿量
                    RO.base.ddx=d_x+ddx;%x方向位置误差矩阵
                    RO.base.ddy=ddy;%y方向位置误差矩阵
                    RO.base.ddz=dz+ddz;%z方向位置误差矩阵
                    RO.base.ddx=d_x+ddx-dx_measure;%x方向位置误差矩阵
                    RO.base.ddz=dz+ddz-dz_measure;%z方向位置误差矩阵 补偿
                
                    RO_err.base.ddx=d_x+ddx;%x方向位置误差矩阵
                    RO_err.base.ddy=ddy;%y方向位置误差矩阵
                    RO_err.base.ddz=dz+ddz;%z方向位置误差矩阵

                    nowfe=my_getDirPtFoc(RO.base);
                    nowfe_err=my_getDirPtFoc(RO_err.base);
                    [maxfe,maxfedir]=max(nowfe);
                    [maxfe_err,maxfedir_err]=max(nowfe_err);
                    RO.feList{loopi,loopj}=RO.feList{loopi,loopj}+nowfe./simLen;%位置、姿态误差
                    dGList(loopi,loopj)=dGList(loopi,loopj)-20*log10(maxfe/max(RO.f));
                    
                    RO_err.feList{loopi,loopj}=RO_err.feList{loopi,loopj}+nowfe_err./simLen;%位置、姿态误差
                    dGList_err(loopi,loopj)=dGList_err(loopi,loopj)-20*log10(maxfe_err/max(RO_err.f));
                    
                    sitaList(loopi,loopj)=sita(maxfedir)*180/pi;%记录当前sita误差
                    sitaList_err(loopi,loopj)=sita(maxfedir_err)*180/pi;
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
            fprintf('%d挠曲指数%f时峰值衰减补偿后：%f dB，补偿前：%f dB,方向角补偿后：%f °，补偿前：%f °\n'...
            ,loopi,2*umax(loopi)*180/pi,dGList(loopi),dGList_err(loopi),sitaList(loopi),sitaList_err(loopi))
        case 2
            fprintf('位置误差标准差：%f，挠曲指数：%f时峰值衰减补偿后：%f dB，补偿前：%f dB,方向角补偿后：%f °，补偿前：%f °\n'...
            ,sigmaList(loopi)/lambda,2*umax(loopi)*180/pi,dGList(loopi,loopi),dGList_err(loopi,loopi),sitaList(loopi,loopi),sitaList_err(loopi,loopi))
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
        plot(umax*180/pi*2,dGList_err,'LineWidth',2);
        hold on;
        plot(umax*180/pi*2,dGList,'LineWidth',2);
        hold on;
        plot(umax((126:375))*180/pi*2,0.5*ones(size(umax((126:375)))),'--','LineWidth',1);
        xlabel('挠曲指数/°','Fontsize',15);
        ylabel('增益损失/dB','Fontsize',15);
        legend("补偿前","补偿后","0.5");
        title('增益损失补偿前后对比图(最大挠曲指数10°,随机安装误差均方差1mm)','Fontsize',15);
        %%
        figure
        
        if dip==5
            sita100=mean(reshape(abs(sitaList(1:700)),100,7),2);
            sita50=sita100(26:75);
            sita50(1:25)=(sita50(1:25)+flipud(sita100(1:25)))/2;
            sita50(26:50)=(sita50(26:50)+flipud(sita100(76:100)))/2;

            sita100_err=mean(reshape(abs(sitaList_err(1:700)),100,7),2);
            sita50_err=sita100_err(26:75);
            sita50_err(1:25)=(sita50_err(1:25)+flipud(sita100_err(1:25)))/2;
            sita50_err(26:50)=(sita50_err(26:50)+flipud(sita100_err(76:100)))/2;
            plot(umax((26:75))*180/pi*2,sita50_err,'LineWidth',2);
            hold on;
            plot(umax((26:75))*180/pi*2,sita50,'LineWidth',2);
        elseif dip ==1
            sita100=mean(reshape(abs(sitaList(1:3500)),500,7),2);
            sita50=sita100(126:375);
            sita50(1:125)=(sita50(1:125)+flipud(sita100(1:125)))/2;
            sita50(126:250)=(sita50(126:250)+flipud(sita100(376:500)))/2;

            sita100_err=mean(reshape(abs(sitaList_err(1:3500)),500,7),2);
            sita50_err=sita100(126:375);
            sita50_err(1:125)=(sita50_err(1:125)+flipud(sita100_err(1:125)))/2;
            sita50_err(126:250)=(sita50_err(126:250)+flipud(sita100_err(376:500)))/2;
            plot(umax((126:375))*180/pi*2,sita50_err,'LineWidth',2);
            hold on;
            plot(umax((126:375))*180/pi*2,sita50,'LineWidth',2);
        end
        
        
        xlabel('挠曲指数/°','Fontsize',15);
        ylabel('方向角误差/°','Fontsize',15);
        legend("补偿前","补偿后");
        title('增益损失补偿前后方向角对比图(最大挠曲指数10°,随机安装误差均方差1mm)','Fontsize',15);
        
        %%
        figure
        
        if dip==5
            sita100=mean(reshape(abs(dGList(1:700)),100,7),2);
            sita50=sita100(26:75);
            sita50(1:25)=(sita50(1:25)+flipud(sita100(1:25)))/2;
            sita50(26:50)=(sita50(26:50)+flipud(sita100(76:100)))/2;

            sita100_err=mean(reshape(abs(dGList_err(1:700)),100,7),2);
            sita50_err=sita100_err(26:75);
            sita50_err(1:25)=(sita50_err(1:25)+flipud(sita100_err(1:25)))/2;
            sita50_err(26:50)=(sita50_err(26:50)+flipud(sita100_err(76:100)))/2;
            plot(umax((26:75))*180/pi*2,sita50_err,'LineWidth',2);
            hold on;
            plot(umax((26:75))*180/pi*2,sita50,'LineWidth',2);
        elseif dip ==1
            sita100=mean(reshape(abs(dGList(1:3500)),500,7),2);
            sita50=sita100(126:375);
            sita50(1:125)=(sita50(1:125)+flipud(sita100(1:125)))/2;
            sita50(126:250)=(sita50(126:250)+flipud(sita100(376:500)))/2;

            sita100_err=mean(reshape(abs(dGList_err(1:3500)),500,7),2);
            sita50_err=sita100(126:375);
            sita50_err(1:125)=(sita50_err(1:125)+flipud(sita100_err(1:125)))/2;
            sita50_err(126:250)=(sita50_err(126:250)+flipud(sita100_err(376:500)))/2;
            plot(umax((126:375))*180/pi*2,sita50_err,'LineWidth',2);
            hold on;
            plot(umax((126:375))*180/pi*2,sita50,'LineWidth',2);
        end
        hold on;
        plot(umax((126:375))*180/pi*2,0.5*ones(size(umax((126:375)))),'--','LineWidth',1);
        xlabel('挠曲指数/°','Fontsize',15);
        ylabel('增益损失/dB','Fontsize',15);
        legend("补偿前","补偿后","0.5");
        title('增益损失补偿前后对比图(最大挠曲指数10°,随机安装误差均方差1mm)','Fontsize',15);
%         title('机翼挠曲变形补偿后增益损失图','Fontsize',15);
%         figure
%         load dGList_noerr5;
%         plot(linspace(0,700,len),dGList_noerr5,'LineWidth',2);
%         hold on;
%         plot(linspace(0,700,len),dGList,'LineWidth',2);
%         xlabel('时间/s','Fontsize',15);
%         ylabel('增益损失/dB','Fontsize',15);
%         legend("补偿前","补偿后");
%         title('增益损失补偿前后对比图(最大挠曲指数5°)','Fontsize',15);
%         xlim([200,500]);
%         title('机翼挠曲变形补偿后增益损失图','Fontsize',15);

%         figure
%         plot(abs(umax(1:len))*5*180/pi,dGList,'LineWidth',2);
%         load dGList_nos;hold on;
%         plot(abs(umax(1:len))*5*180/pi,dGList_nos,'LineWidth',2);
%         xlabel('挠曲指数/°','Fontsize',15);
%         ylabel('增益损失/dB','Fontsize',15);
%         legend("补偿前","补偿后");
%         title('增益损失补偿前后对比图(最大挠曲指数5°)','Fontsize',15);
%         title('机翼挠曲变形补偿后增益损失图','Fontsize',15);
    case 2
        figure
        x=repmat(sigmaList/lambda,len,1);
        y=repmat(umax'*180/pi,1,len);
        z=dGList';
        surf(x,y,z);
        xlabel('位置随机误差均方差/\lambda','Fontsize',15);
        ylabel('挠曲角/°','Fontsize',15);
        zlabel('增益损失/dB','Fontsize',15);
        title('机翼挠曲变形+随机误差增益损失图','Fontsize',15);
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
