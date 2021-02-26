%% �����λ����������ֵ˥�� 14��������޸�
% �汾ʱ�䣺2020/01/21 
% ��Ƴ��ԣ��������������ζ���Ƥ���ߵ�Ӱ�� ��װ 2m�����İ�װ��Ƥ ��ϸ�ؽ���λ��
% ����ͼ��flag 1 200 s - 500 s �Ĳ������ - ���ֵ -��ֵ
% ����ͼ��flag 2 200 s - 500 s �Ĳ������ - ÿһ������ʱ�̶�һ�������װ���ǽ��ж�����ؿ�������������ؿ���ľ�ֵ
% ���� 200 s - 500 s �����ֵ ��ֵ
clc
clear
close all
sita=linspace(-pi/2,pi/2,2000);
phi=0;
cm=0.01;
M=33;%x������Ԫ����
N=21;%y��������M
lambda=5*cm;%����
dx=lambda*0.6;%x������Ԫ����
dy=lambda*0.6;%y������Ԫ����
xm=floor(M/2)*dx;%x����������
ym=floor(N/2)*dy;%y����������
RO.base.dx=lambda*0.6;%x������Ԫ����
RO.base.dy=lambda*0.6;%y������Ԫ����
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
simLen=1;%���ؿ���������
Zmax=1*(0:0.025:0.5)*lambda;%���Ʊ���ָ��
% umax=Zmax/max(Zmax)*5*pi/180;%���Ʊ���ָ��
ts=0.1;%����ʱ��
nts=2*ts;
samp.begin=200/nts;
samp.end=500/nts;
samp.dir=10*2;
load NSINSFR3;
% umax=NSINSFR.attnewall_true(2,1:2:7000)/2;%��λ������������ʵֵ �غ�ģ��
% umax_measure=NSINSFR.attnewall(2,1:2:7000)/2;%��λ��������������ֵ �غ�ģ��
% umax=NSINSFR.attnewall_true(2,:)/2;%��λ������������ʵֵ �غ�ģ��
% umax_measure=NSINSFR.attnewall(2,:)/2;%��λ��������������ֵ �غ�ģ��
umax=NSINSFR.attnewall_true(2,samp.begin:samp.dir:samp.end)/2;%��λ������������ʵֵ �غ�ģ��
umax_measure=NSINSFR.attnewall(2,samp.begin:samp.dir:samp.end)/2;%��λ��������������ֵ �غ�ģ��
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
        fprintf('����...%6.1f %%%6.1f %%',0,0);
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
                xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%��Ƥ���İ�װ�� 2m ��
                [d_x,dz]=my_fixdefload(xdir+ddx,1,umax(loopi),1);%�����غ�����ģ�ͱ�׼�� ��׼ģ��1
                [dx_measure,dz_measure]=my_fixdefload(xdir+ddx,1,umax_measure(loopi),1);%�����غ�����ģ�Ͳ�����
                RO.base.ddx=d_x;%x����λ��������
                RO.base.ddy=ddy;%y����λ��������
                RO.base.ddz=dz+ddz;%z����λ��������
                RO.base.ddx=d_x-dx_measure;%x����λ��������
                RO.base.ddz=dz+ddz-dz_measure;%z����λ�������� ����
                
                nowfe=my_getDirPtFoc(RO.base);
                RO.feList{loopi}=RO.feList{loopi}+nowfe./simLen;%λ�����
                dGList(loopi)=dGList(loopi)-20*log10(max(nowfe)/max(RO.f));
            case 2
                for loopj=1:len_loopj
                    ddx=sigmaList(loopi)*randn(N,M);
                    ddy=sigmaList(loopi)*randn(N,M);
                    ddz=sigmaList(loopi)*randn(N,M);
                    xdir=2 + (-floor(M/2):1:floor(M/2))*dx;%��Ƥ���İ�װ�� 2m ��
                    [d_x,dz]=my_fixdefload(xdir+ddx,1,umax(loopj),1);%�����غ�����ģ�ͱ�׼�� ��׼ģ��1
                    [dx_measure,dz_measure]=my_fixdefload(xdir+ddx,1,umax_measure(loopj),1);%�����غ�����ģ�Ͳ�����
                    RO.base.ddx=d_x;%x����λ��������
                    RO.base.ddy=ddy;%y����λ��������
                    RO.base.ddz=dz+ddz;%z����λ��������
                    RO.base.ddx=d_x-dx_measure;%x����λ��������
                    RO.base.ddz=dz+ddz-dz_measure;%z����λ�������� ����
                
%                 RO.base.ddz=dz+ddz-dz_measure;%z����λ�������� ����
                
                    nowfe=my_getDirPtFoc(RO.base);
                    RO.feList{loopi,loopj}=RO.feList{loopi,loopj}+nowfe./simLen;%λ�á���̬���
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
            fprintf('%d����ָ��%fʱ��ֵ˥����%f dB\n'...
            ,loopi,5*umax(loopi)*180/pi,dGList(loopi))
        case 2
            fprintf('λ������׼�%f����ֵ˥�����ֵ��%f dB ��ֵ��%f dB\n'...
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
legend('�����','������');
xlabel('theta/radian');
ylabel('amplitude');
title('��һ������ͼphi=0')
%% ��ά����ͼ

%%
switch flag
    case 1
        %��ͼ
        figure
        plot(linspace(samp.begin*nts,samp.end*nts,len_loopi),dGList,'LineWidth',2);
        xlabel('ʱ��/s','Fontsize',15);
        ylabel('������ʧ/dB','Fontsize',15);
%         title('������������δ����������ʧͼ','Fontsize',15);
%         title('�����������β�����������ʧͼ(\sigma=0.001mm)','Fontsize',15);
        title('�����������β�����������ʧͼ','Fontsize',15);
        fprintf("���������ʧ%fdB��ƽ��������ʧ%fdB\n",max(dGList),mean(dGList));
        %��ͼ
        figure
        plot(linspace(samp.begin*nts,samp.end*nts,len_loopi),dGList,'LineWidth',2);
        hold on;
        plot(linspace(samp.begin*nts,samp.end*nts,len_loopi),max(dGList)*ones(1,len_loopi),"-.",'LineWidth',2);
        hold on;
        plot(linspace(samp.begin*nts,samp.end*nts,len_loopi),mean(dGList)*ones(1,len_loopi),"--",'LineWidth',2);
        xlabel('ʱ��/s','Fontsize',15);
        ylabel('������ʧ/dB','Fontsize',15);
        title('������������δ����������ʧͼ','Fontsize',15);
%         title('�����������β�����������ʧͼ(\sigma=0.001mm)','Fontsize',15);
        legend("������","���ֵ","��ֵ");
    case 2
        figure
        x=repmat(sigmaList/lambda,len_loopj,1);
        y=repmat(linspace(samp.begin*nts,samp.end*nts,len_loopj)',1,len_loopi);
        z=dGList';
        surf(x,y,z);
        xlabel('λ�������������/\lambda','Fontsize',15);
        ylabel('ʱ��/s','Fontsize',15);
        zlabel('������ʧ/dB','Fontsize',15);
        title('������������+������������ʧͼ','Fontsize',15); 
        figure
        plot(sigmaList/lambda,simResult(:,1),'LineWidth',2);
        hold on;
        plot(sigmaList/lambda,simResult(:,2),'LineWidth',2);
        xlabel('��װλ����������/\lambda','Fontsize',15);
        ylabel('������ʧ/dB','Fontsize',15);
        title('��λ��Ч��','Fontsize',15);
        legend("���ֵ","ƽ��ֵ");
end
%%
figure
ydir=(-floor(N/2):1:floor(N/2))*dy;
[MList,NList]=meshgrid(xdir,ydir);
[d_x,dz]=my_fixdefload(xdir+ddx,5,1*pi/180,1);%�����غ�����ģ�ͱ�׼�� ��׼ģ��1
switch mod(flag,3)
    case 1
        subplot(211)
        surf((MList+d_x+0*ddx)/cm,(NList+0*ddy)/cm,(dz+0*ddz)/cm);
%         title('������������ͼ','Fontsize',15);
        title('����ָ��1��ʱ������Ԫλ��','Fontsize',15);
        ylim([-50,50]);zlim([0,10]);
        colorbar;
        subplot(212)
        surf((MList+d_x+ddx)/cm,(NList+0*ddy)/cm,(dz+ddz)/cm);
%         title('������������ͼ','Fontsize',15);
        title('����ָ��1��ʱ������Ԫλ��(\sigma=0.001mm)','Fontsize',15);
        ylim([-50,50]);zlim([0,10]);
        colorbar;
end
xlabel('x/cm','Fontsize',15);ylabel('y/cm','Fontsize',15);zlabel('z/cm','Fontsize',15);
