%% �����λ����������ֵ˥��
% �汾ʱ�䣺2020/12/20 
% ��Ƴ��ԣ��������������ζ���Ƥ���ߵ�Ӱ�� ��װ 2m�����İ�װ��Ƥ �����ؽ���λ��
clc
clear
close all
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

RO.base.phi=phi;
RO.base.sita=sita;
RO.base.lambda=lambda;%����
RO.base.M=M;%x������Ԫ����
RO.base.N=N;%y������Ԫ����
RO.base.sita0=0;%��������ָ������
RO.base.phi0=0;%��������ָ��λ��
RO.base.ddsita=0;%������������
RO.base.ddphi=0;%��λ��������
RO.base.ddx=zeros(N,M);%x����λ��������
RO.base.ddy=zeros(N,M);%y����λ��������
RO.base.ddz=zeros(N,M);%z����λ��������

RO.f=my_getDirPtFoc(RO.base);%��λ�����
flag=1;%���Ʒ���ģʽ 1��x������������� 2����������+�����άλ����� 
SigmaMax=0.1*lambda;%���Ʒ���λ��������ֵ
SitaMax=1;%���Ʒ��温��������ֵ
PhiMax=1;%���Ʒ��淽��������ֵ
simLen=100;%���ؿ���������
Zmax=1*(0:0.025:0.5)*lambda;%���Ʊ���ָ��
umax=Zmax/max(Zmax)*0.5*pi/180;%���Ʊ���ָ��
sigmaList=Zmax/max(Zmax)*SigmaMax;
sigmasitaList=Zmax/max(Zmax)*SitaMax*pi/180;
sigmaphiList=Zmax/max(Zmax)*PhiMax*pi/180;
if flag<2
    dGList=zeros(1,length(Zmax));
    RO.feList=cell(size(Zmax));
    for loopi=1:length(Zmax)
        RO.feList{loopi}=zeros(size(RO.f));
    end
else
    dGList=zeros(length(Zmax),length(sigmasitaList));
    RO.feList=cell(length(Zmax),length(sigmasitaList));
    for loopi=1:length(Zmax)
        for loopj=1:length(sigmasitaList)
            RO.feList{loopi,loopj}=zeros(size(RO.f));
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
                xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%��Ƥ���İ�װ�� 2m ��
                dz=my_fixdefload(xdir,5,umax(loopi));%�����غ�����ģ�ͱ�׼��
                dz_measure=my_fixdefload(xdir,5,umax(loopi)-0.1*pi/180);%�����غ�����ģ�Ͳ�����
                RO.base.ddx=ddx;%x����λ��������
                RO.base.ddy=ddy;%y����λ��������
                RO.base.ddz=dz+ddz;%z����λ��������
                
%                 RO.base.ddz=dz+ddz-dz_measure;%z����λ�������� ����
                
                nowfe=my_getDirPtFoc(RO.base);
                RO.feList{loopi}=RO.feList{loopi}+nowfe./simLen;%λ�����
                dGList(loopi)=dGList(loopi)-20*log10(max(nowfe)/max(RO.f));
            case 2
                for loopj=1:length(Zmax)
                    ddx=sigmaList(loopi)*randn(N,M);
                    ddy=sigmaList(loopi)*randn(N,M);
                    ddz=sigmaList(loopi)*randn(N,M);
                    xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%��Ƥ���İ�װ�� 2m ��
                    dz=my_fixdefload(xdir,5,umax(loopj));%�����غ�����ģ�ͱ�׼��
                    dz_measure=my_fixdefload(xdir,5,umax(loopj)-0.1*pi/180);%�����غ�����ģ�Ͳ�����
                    RO.base.ddx=ddx;%x����λ��������
                    RO.base.ddy=ddy;%y����λ��������
                    RO.base.ddz=dz+ddz;%z����λ��������
                
%                 RO.base.ddz=dz+ddz-dz_measure;%z����λ�������� ����
                
                    nowfe=my_getDirPtFoc(RO.base);
                    RO.feList{loopi,loopj}=RO.feList{loopi,loopj}+nowfe./simLen;%λ�á���̬���
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
        for loopj=1:length(Zmax)
            dGList(loopi,loopj)=dGList(loopi,loopj)/simLen;
        end
    end
    fprintf(' ');
    switch flag
        case 1 
            fprintf('������%fʱ��ֵ˥����%f dB\n'...
            ,umax(loopi)*180/pi,dGList(loopi))
        case 2
            fprintf('λ������׼�%f�������ǣ�%fʱ��ֵ˥����%f dB\n'...
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
legend('�����','������');
xlabel('theta/radian');
ylabel('amplitude');
title('��һ������ͼphi=0')
%% ��ά����ͼ

%%
switch flag
    case 1
        %���
        [xi,yi,fuc]=my_fit(umax*180/pi,dGList);
        %��ͼ
        figure
        plot(umax*180/pi,dGList,'*');
        xlabel('������/��','Fontsize',15);
        hold on;plot(xi,yi);
        ylabel('������ʧ/dB','Fontsize',15);
        title('������������������ʧͼ','Fontsize',15);
        legend('����','���');
    case 2
        figure
        x=repmat(sigmaList/lambda,length(Zmax),1);
        y=repmat(umax'*180/pi,1,length(sigmaList));
        z=dGList';
        surf(x,y,z);
        xlabel('λ�������������/\lambda','Fontsize',15);
        ylabel('������/��','Fontsize',15);
        zlabel('������ʧ/dB','Fontsize',15);
        title('������������+������������ʧͼ','Fontsize',15);
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
