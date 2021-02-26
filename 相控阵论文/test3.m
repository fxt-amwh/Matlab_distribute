%% �����λ����������ֵ˥��
% �汾ʱ�䣺2020/12/19 
% ��Ƴ��ԣ��������������ζ���Ƥ���ߵ�Ӱ��
clc
clear
close all
global lambda
global M
global N
sita=linspace(-pi/2,pi/2,2000);
phi=-pi/2;
cm=0.01;
M=33;%x������Ԫ����
N=9;%y��������M
lambda=3*cm;%����
dx=lambda*2;%x������Ԫ����
dy=lambda*2;%y������Ԫ����
xm=floor(M/2)*dx;%x����������
ym=floor(N/2)*dy;%y����������
f=my_getDirPt(sita,phi,0,0,0);%��λ�����
flag=1;%���Ʒ���ģʽ 1��x������������� 2����������+�����άλ����� 
SigmaMax=0.1*lambda;%���Ʒ���λ��������ֵ
SitaMax=1;%���Ʒ��温��������ֵ
PhiMax=1;%���Ʒ��淽��������ֵ
simLen=100;%���ؿ���������
Zmax=1*(0:0.025:0.5)*lambda;%���Ʊ���ָ��
% Zmax=1*(0:0.02:0.4)*lambda;%���Ʊ���ָ��
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
         fprintf('������');
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
                feList{loopi}=feList{loopi}+nowfe./simLen;%λ�����
                dGList(loopi)=dGList(loopi)-20*log10(max(nowfe)/max(f));
            case 2
                for loopj=1:length(Zmax)
                    ddx=sigmaList(loopi)*randn(N,M);
                    ddy=sigmaList(loopi)*randn(N,M);
                    ddz=sigmaList(loopi)*randn(N,M);
                    xdir=(-floor(M/2):1:floor(M/2))*dx;
                    dz=my_fixdef(xdir,xm,ym,Zmax(loopj),1,ddx,ddy);
                    nowfe=my_getDirPt(sita,phi,ddx,ddy,dz+ddz);
                    feList{loopi,loopj}=feList{loopi,loopj}+nowfe./simLen;%λ�á���̬���
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
            fprintf('����ָ��%fʱ��ֵ˥����%f dB\n'...
            ,Zmax(loopi)/lambda,dGList(loopi))
        case 2
            fprintf('λ������׼�%f������ָ����%fʱ��ֵ˥����%f dB\n'...
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
legend('�����','������');
xlabel('theta/radian');
ylabel('amplitude');
title('��һ������ͼphi=0')
%% ��ά����ͼ

%%
switch flag
    case 1
        %���
        [xi,yi,fuc]=my_fit(Zmax/lambda,dGList);
        %��ͼ
        figure
        plot(Zmax/lambda,dGList,'*');
        xlabel('����ָ��/\lambda','Fontsize',15);
        hold on;plot(xi,yi);
        ylabel('������ʧ/dB','Fontsize',15);
        title('��������������ʧͼ','Fontsize',15);
        legend('����','���');
    case 2
        figure
        x=repmat(sigmaList/lambda,length(Zmax),1);
        y=repmat(Zmax'/lambda,1,length(sigmaList));
        z=dGList';
        surf(x,y,z);
        xlabel('λ�������������/\lambda','Fontsize',15);
        ylabel('����ָ��/\lambda','Fontsize',15);
        zlabel('������ʧ/dB','Fontsize',15);
        title('��������+������������ʧͼ','Fontsize',15);
end
if flag<2
    disp('��Ϻ���:')
    fuc
end
%%
figure
ydir=(-floor(N/2):1:floor(N/2))*dy;
[MList,NList]=meshgrid(xdir,ydir);
switch mod(flag,3)
    case 1
        surf((MList+ddx)/cm,(NList+ddy)/cm,(dz+ddz)/cm);
        title('��������','Fontsize',15);
end
xlabel('x/cm','Fontsize',15);ylabel('y/cm','Fontsize',15);zlabel('z/cm','Fontsize',15);
