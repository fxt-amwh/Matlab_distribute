%% һ�� + һ�� ��Ҫ
%�ļ����ƣ�test10
%������Ϣ��������
%������������Ե����˲� �ٶ�ƥ��� �з��� ���ؿ������
%�汾ʱ�䣺2021/3/15 21:22
% 1��������ƹ켣�˶����������ߴ��������
% 2��������ƻ����˶��������ӹߴ���������
% 3���ӹ�����ע�����
% 4����������Ե���
% 5���˲�
% 6��������
%% �����˶�������Ʒ���
clc
clear
close all
simLen=100;
gvar;    % ����ȫ�ֱ���
ts = 0.1;%����ʱ��
att0 = [0;0;0]*arcdeg; vn0 = [0;100;0]; pos0 = [[34;108]*arcdeg;100];
%     ���������� ��������� ��λ������ ������ٶ� ����ʱ�� �����ι켣��ƣ�
wat = [  0,         0,          0,         0,         10       %��ֹ
        60,         0,          0,         0,         60       %����
        -60,        0,          0,         0,         60       %����
        0,          0,          0,         0,         50       %����
        0,          0,          60,       0,         360      %����ת��
        0,          0,          0,         0,         50       %����
        -60,        0,          0,         0,         60       %����
        60,         0,          0,         0,        60       %����
        0,          0,          0,         0,         10       %����
        ];    %��ֹ
wat(:,1:3) = wat(:,1:3)*arcdeg/60;  % deg/min->deg/s
fprintf('�켣������...'); 
[attm, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
fprintf('�켣������ɣ�\n���߷�����...');
[wm, vm] = my_av2imu(attm, vn, pos, ts);
tt = (0:length(attm)-1)'*ts;
%% ������Ϣ
wim0=my_getII(wm(1,:)',ts);
wim=[wim0,my_getII(wm',ts)];%wimע������������ϵ�µ�����
fm=[my_getII(vm(1,:)',ts),my_getII(vm',ts)];
MINS.wim=wim;
MINS.wim0=wim0;
MINS.fm=fm;
fprintf('���߷�����ɣ�\n');
%% �ӹ��˶��������
fprintf('�ӹ߷���...'); 
len=size(wm,1)+1;% len=7200;
um=0*arcdeg;
f=0.01;
u=[zeros(1,len);um*sin(2*pi*f*(0:ts:((len-1)*ts)));zeros(1,len)];
R0=[2;0;0];%����ʱ ���Ӽ�0ʱ��ǰ��ʼ���λ��
Rf0=[0;0;0];%����ʱ ��0ʱ��ǰƫ��
Rf=[-R0(1)*sin(u(2,:)).*sin(u(2,:));zeros(1,len);R0(1)*sin(u(2,:)).*cos(u(2,:))];
R=R0+Rf;
att=2*u;
att_s0=[0;0;0]*arcdeg;%����ʱ ���Ӽ�0ʱ��ǰ��ʼ��ԽǶ�
%% �ӹ߷���������
atterr0=0*[1;2;3]*arcdeg;%ע��������ת˳��zxy
% atterr0=[0;0;0]*arcdeg;%ע��������ת˳��zxy
Cmserr=my_a2mat(atterr0);%ע���м䴦����ת˳��yxz
diffRf0=[0;0;0];%Rf�ĳ�ʼ�仯��
U0=diffRf0+cross(wim0,R0);
SINS1=my_invRI(R,att,atterr0,wim,fm,ts,R0,att_s0,U0);
fprintf('�ӹ߷�����ɣ�\n���ӹ���Ե���...%5.0f %%',0);  
%% �ӹ�������
SERR1.eb = [1;1;1]*dph; SERR1.web = [1;1;1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR1.db = [200;200;200]*ug; SERR1.wdb = [200;200;200]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
% SERR1.eb = 0*[1;1;1]*dph; SERR1.web = 0*[1;1;1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
% SERR1.db = 0*[200;200;200]*ug; SERR1.wdb = 0*[200;200;200]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������

%%
errall1=zeros(3,len/2);
errall2=zeros(3,len/2);
RMSE1=zeros(4,simLen);
RMSE2=zeros(4,simLen);
for loopi=1:simLen
    %% ��Ե���
    atterr0_1=0*[2;3;2]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
    SINSR1=my_SINSgetResult(MINS,SINS1,len/2,SERR1,atterr0_1);
%     fprintf('������Ե�����ɣ�\n�˲�...%5.0f %%',0);
    %% �˲�
    nts = 2*ts; %�������Ͳ���ʱ��
    KFinit1.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
    KFinit1.rk = [0.001;0.001;0.001];  
    KFinit1.Rk = diag(KFinit1.rk)^2;
    KFinit1.P0 = diag([[0.1;0.1;0.1]*arcdeg; [0.001;0.001;0.001]; [0.001;0.001;0.001];
             [1;1;1]*dph; [200;200;200]*ug])^2;
    [SINSFR1,~]=my_getFResult(MINS,SINS1,KFinit1,atterr0_1,len/2,...
        SINSR1.ws_m_addnoise,SINSR1.fs_m_addnoise);
    fprintf('�˲���ɣ�\n'); 
    %% �Ƕ���ֵ zxyת��
    atttrue=zeros(3,len);
    for i=1:len%��ʵ�Ƕ�
       [rz,rx,ry]=dcm2angle(SINS1.Cms{i},'zxy');
       atttrue(:,i)=[rx,ry,-rz]; 
    end
    %% �ٶ�ƥ�� ������
    MERR.eb=[0;0;0]*dph;
    MERR.web=[0;0;0]*dpsh;
    MERR.db=[0;0;0]*ug;
    MERR.wdb=[0;0;0]*ugpsHz;
    Flag.EnBack=true;
    SINS_Ret_VLF=zeros(floor(len/2),22);kk=1;t=0;
    qnb_vm=a2qua(att0);
    vn_vm=vn0;
    posMINS_vm=pos0;
    qnb_vs=a2qua(att0+atterr0_1);
    vn_vs=vn0;
    posSINS_vs=pos0;
    cl = cos(posSINS_vs(1,1)); Re = 6378137;
    posSINS_vs(2,1)=posSINS_vs(2,1)-R0(1)/(cl*Re);
    Filter.X=zeros(12,len/2);
    Cmserr_1=a2mat(att0+atterr0_1);%ע���м䴦����ת˳��yxz ��׼ȷ�ĳ�ʼ��
    qns=m2qua(Cmserr_1);%��ʼ��̬ ��׼ȷ�ĳ�ʼ��
    fs0=[0;0;0];
    KFinit.Qk = diag([SERR1.wdb; SERR1.web;])^2*nts;
    KFinit.rk = [0.001;0.001;0.001];  
    KFinit.Rk = diag(KFinit.rk)^2;
    KFinit.P0 = diag([[0.001;0.001;0.001];[0.1;0.1;0.1]*arcdeg;SERR1.db; SERR1.eb])^2;
    eth0 = earth(posMINS_vm, vn_vm);  % ������ز�������
    kfft=my_kfftV12(eth0,q2mat(qns),fs0,nts);
    kf = kfinit(KFinit.Qk, KFinit.Rk, KFinit.P0,kfft.phi,kfft.H);  % kf�˲�����ʼ��
    % SERR.eb=[0;0;0]*dph;
    % SERR.web=[0;0;0]*dpsh;
    % SERR.db=[0;0;0]*ug;
    % SERR.wdb=[0;0;0]*ugpsHz;
    for k=1:len/2
        t=t+nts;
        [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
            MINS.fm(:,(2*k-1):(2*k))'*ts, ...
            MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
        [ws1,vs1]=imuadderr(SINS1.wis(:,(2*k-1):(2*k))'*ts, ...
            SINS1.fs(:,(2*k-1):(2*k))'*ts, ...
            SERR1.eb, SERR1.web, SERR1.db, SERR1.wdb, ts);
        [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL,ethm]=...
            my_Vrelinsupdate...
            (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R0,ts);

        fs=sum(SINS1.fs(:,(2*k-1):(2*k)),2)/2;
        kfft=my_kfftV12(ethm,q2mat(qnb_vs),fs,nts);
        kf.Phikk_1=kfft.phi;
        kf.Gammak=kfft.Gammak;

        ZV=vn_vs-(vn_vm+VL);% ��������

        kf = kfupdate(kf,ZV,'B');%�������˲�

        if Flag.EnBack
            qnb_vs = qdelphi(qnb_vs,kf.Xk(4:6));  kf.Xk(4:6) = 0;  % ������
            vn_vs = vn_vs - kf.Xk(1:3);  kf.Xk(1:3) = 0;
        end

        Filter.X(:,k)=kf.Xk;
        Filter.X(4:6,k)=q2att(rv2q(Filter.X(4:6,k)));
        SINS_Ret_VLF(kk,:)=...
            [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
        kk=kk+1;
%         if mod(t,100)<nts,disp(fix(t));end
    end
    Filter.XT=Filter.X';
    
    %% ͳ��RMSE
    beg_T=200;
    end_T=500;
    beg_index=beg_T/nts;
    end_index=end_T/nts;
    %diag(resultR(:,sam.begin:sam.end)*resultR(:,sam.begin:sam.end)')/length(resultR(:,sam.begin:sam.end))/simLen^2
    errall1=errall1+abs(SINSFR1.attall(:,1:end)-atttrue(:,1:2:end));
    err1=errall1(:,beg_index:end_index)/loopi;
    RMSE1(1:3,loopi) = 1000*sqrt(diag(err1*err1'/length(err1)));
  
    if  Flag.EnBack
        tmp=abs(SINS_Ret_VLF(1:end,1:3)'-SINS_Ret_VLF(1:end,10:12)'+[atttrue(1:2,1:2:end);-atttrue(3,1:2:end)]);
        tmp(3,tmp(3,:)<-pi)=0.2*pi/180;tmp(3,tmp(3,:)>pi)=0.2*pi/180;
        errall2=errall2+abs(tmp);
        err2=errall2(:,beg_index:end_index)/loopi;
        err2(3,err2(3,:)<-pi)=0.1*pi/180;err2(3,err2(3,:)>pi)=0.1*pi/180;
        RMSE2(1:3,loopi) = 1000*sqrt(diag(err2*err2'/length(err2)));
    else    
        errall2=errall2+abs(Filter.XT(:,4:6)+atterr0_1')';
        err2=errall2(:,beg_index:end_index)/loopi;
        RMSE2(1:3,loopi) = 1000*sqrt(diag(err2*err2'/length(err2)));
    end
    RMSE1(4,loopi)=sqrt(RMSE1(1:3,loopi)'*RMSE1(1:3,loopi));
    RMSE2(4,loopi)=sqrt(RMSE2(1:3,loopi)'*RMSE2(1:3,loopi));
    fprintf("loop: %d RMSE��Ե���:(mrad):%f %f %f %f ",loopi,RMSE1(:,loopi)')
    fprintf("RMSE�ٶ�ƥ��:(mrad)::%f %f %f %f \r\n",RMSE2(:,loopi)');
end
%% ͳ��RMSE
fprintf("RMSE��Ե���:(mrad)");
disp(RMSE1(:,simLen)');
fprintf("RMSE�ٶ�ƥ��:(mrad)");
disp(RMSE2(:,simLen)');
% Mxyz=my_Lhtoxyz_all(SINS_Ret_VLF(:,7:9));
% Sxyz=my_Lhtoxyz_all(SINS_Ret_VLF(:,16:18));
% difxyz=Sxyz-Mxyz;
% plot(diag(sqrt(difxyz*difxyz')),'DisplayName','difxyz')