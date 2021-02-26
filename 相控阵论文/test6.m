%% �����λ����������ֵ˥��
% �汾ʱ�䣺2021/1/23 
% ��Ƴ��ԣ��������������ζ���Ƥ���ߵ�Ӱ�� ��װ 2m�����İ�װ��Ƥ ��ϸ�ؽ���λ��
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
RO.base.dx=lambda*0.6;%x������Ԫ����
RO.base.dy=lambda*0.6;%y������Ԫ����
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
% umax=Zmax/max(Zmax)*5*pi/180;%���Ʊ���ָ��
load NSINSFR8;
% umax=NSINSFR.attnewall_true(2,1:2:7000)/2;%��λ������������ʵֵ �غ�ģ��
% umax_measure=NSINSFR.attnewall(2,1:2:7000)/2;%��λ��������������ֵ �غ�ģ��
dip=1;
umax=NSINSFR.attnewall_true(2,1:dip:end)/2;%��λ������������ʵֵ �غ�ģ��
umax_measure=NSINSFR.attnewall(2,1:dip:end)/2;%��λ��������������ֵ �غ�ģ��

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
         fprintf('������');
    end
    for simi=1:simLen
        switch flag
            case 1
                ddx=0.001*randn(N,M);
                ddy=0.001*randn(N,M);
                ddz=0.001*randn(N,M);
                xdir=0.5 + (-floor(M/2):1:floor(M/2))*dx;%��Ƥ���İ�װ�� 2m ��
                [d_x,dz]=my_fixdefload(xdir,2,2*umax(loopi),1);%�����غ�����ģ�ͱ�׼�� ��׼ģ��1
                [dx_measure,dz_measure]=my_fixdefload(xdir,2,2*umax_measure(loopi),1);%�����غ�����ģ�Ͳ�����
                
                RO.base.ddx=d_x+ddx-dx_measure;%x����λ��������
                RO.base.ddy=ddy;%y����λ��������
                RO.base.ddz=dz+ddz-dz_measure;%z����λ�������� ����
                
                RO_err.base.ddx=d_x+ddx;%x����λ��������
                RO_err.base.ddy=ddy;%y����λ��������
                RO_err.base.ddz=dz+ddz;%z����λ��������
                
                nowfe=my_getDirPtFoc(RO.base);%������
                nowfe_err=my_getDirPtFoc(RO_err.base);%������
                RO.feList{loopi}=RO.feList{loopi}+nowfe./simLen;%λ�����
                RO_err.feList{loopi}=RO_err.feList{loopi}+nowfe_err./simLen;%λ�����
                [maxfe,maxfedir]=max(nowfe);
                [maxfe_err,maxfedir_err]=max(nowfe_err);
                dGList(loopi)=dGList(loopi)-20*log10(maxfe/max(RO.f));
                dGList_err(loopi)=dGList_err(loopi)-20*log10(maxfe_err/max(RO_err.f));
                
                sitaList(loopi)=sita(maxfedir)*180/pi;%��¼��ǰsita���
                sitaList_err(loopi)=sita(maxfedir_err)*180/pi;
            case 2
                for loopj=1:len
                    ddx=sigmaList(loopi)*randn(N,M);
                    ddy=sigmaList(loopi)*randn(N,M);
                    ddz=sigmaList(loopi)*randn(N,M);
                    xdir=0.5 + (-floor(M/2):1:floor(M/2))*dx;%��Ƥ���İ�װ�� 2m ��
                    [d_x,dz]=my_fixdefload(xdir,2,2*umax(loopj),1);%�����غ�����ģ�ͱ�׼�� ��׼ģ��1
                    [dx_measure,dz_measure]=my_fixdefload(xdir,2,2*umax_measure(loopj),1);%�����غ�����ģ�Ͳ�����
                    RO.base.ddx=d_x+ddx;%x����λ��������
                    RO.base.ddy=ddy;%y����λ��������
                    RO.base.ddz=dz+ddz;%z����λ��������
                    RO.base.ddx=d_x+ddx-dx_measure;%x����λ��������
                    RO.base.ddz=dz+ddz-dz_measure;%z����λ�������� ����
                
                    RO_err.base.ddx=d_x+ddx;%x����λ��������
                    RO_err.base.ddy=ddy;%y����λ��������
                    RO_err.base.ddz=dz+ddz;%z����λ��������

                    nowfe=my_getDirPtFoc(RO.base);
                    nowfe_err=my_getDirPtFoc(RO_err.base);
                    [maxfe,maxfedir]=max(nowfe);
                    [maxfe_err,maxfedir_err]=max(nowfe_err);
                    RO.feList{loopi,loopj}=RO.feList{loopi,loopj}+nowfe./simLen;%λ�á���̬���
                    dGList(loopi,loopj)=dGList(loopi,loopj)-20*log10(maxfe/max(RO.f));
                    
                    RO_err.feList{loopi,loopj}=RO_err.feList{loopi,loopj}+nowfe_err./simLen;%λ�á���̬���
                    dGList_err(loopi,loopj)=dGList_err(loopi,loopj)-20*log10(maxfe_err/max(RO_err.f));
                    
                    sitaList(loopi,loopj)=sita(maxfedir)*180/pi;%��¼��ǰsita���
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
            fprintf('%d����ָ��%fʱ��ֵ˥��������%f dB������ǰ��%f dB,����ǲ�����%f �㣬����ǰ��%f ��\n'...
            ,loopi,2*umax(loopi)*180/pi,dGList(loopi),dGList_err(loopi),sitaList(loopi),sitaList_err(loopi))
        case 2
            fprintf('λ������׼�%f������ָ����%fʱ��ֵ˥��������%f dB������ǰ��%f dB,����ǲ�����%f �㣬����ǰ��%f ��\n'...
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
        plot(umax*180/pi*2,dGList_err,'LineWidth',2);
        hold on;
        plot(umax*180/pi*2,dGList,'LineWidth',2);
        hold on;
        plot(umax((126:375))*180/pi*2,0.5*ones(size(umax((126:375)))),'--','LineWidth',1);
        xlabel('����ָ��/��','Fontsize',15);
        ylabel('������ʧ/dB','Fontsize',15);
        legend("����ǰ","������","0.5");
        title('������ʧ����ǰ��Ա�ͼ(�������ָ��10��,�����װ��������1mm)','Fontsize',15);
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
        
        
        xlabel('����ָ��/��','Fontsize',15);
        ylabel('��������/��','Fontsize',15);
        legend("����ǰ","������");
        title('������ʧ����ǰ����ǶԱ�ͼ(�������ָ��10��,�����װ��������1mm)','Fontsize',15);
        
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
        xlabel('����ָ��/��','Fontsize',15);
        ylabel('������ʧ/dB','Fontsize',15);
        legend("����ǰ","������","0.5");
        title('������ʧ����ǰ��Ա�ͼ(�������ָ��10��,�����װ��������1mm)','Fontsize',15);
%         title('�����������β�����������ʧͼ','Fontsize',15);
%         figure
%         load dGList_noerr5;
%         plot(linspace(0,700,len),dGList_noerr5,'LineWidth',2);
%         hold on;
%         plot(linspace(0,700,len),dGList,'LineWidth',2);
%         xlabel('ʱ��/s','Fontsize',15);
%         ylabel('������ʧ/dB','Fontsize',15);
%         legend("����ǰ","������");
%         title('������ʧ����ǰ��Ա�ͼ(�������ָ��5��)','Fontsize',15);
%         xlim([200,500]);
%         title('�����������β�����������ʧͼ','Fontsize',15);

%         figure
%         plot(abs(umax(1:len))*5*180/pi,dGList,'LineWidth',2);
%         load dGList_nos;hold on;
%         plot(abs(umax(1:len))*5*180/pi,dGList_nos,'LineWidth',2);
%         xlabel('����ָ��/��','Fontsize',15);
%         ylabel('������ʧ/dB','Fontsize',15);
%         legend("����ǰ","������");
%         title('������ʧ����ǰ��Ա�ͼ(�������ָ��5��)','Fontsize',15);
%         title('�����������β�����������ʧͼ','Fontsize',15);
    case 2
        figure
        x=repmat(sigmaList/lambda,len,1);
        y=repmat(umax'*180/pi,1,len);
        z=dGList';
        surf(x,y,z);
        xlabel('λ�������������/\lambda','Fontsize',15);
        ylabel('������/��','Fontsize',15);
        zlabel('������ʧ/dB','Fontsize',15);
        title('������������+������������ʧͼ','Fontsize',15);
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
