%�ļ����ƣ�test2
%������Ϣ��������
%�����������ٶ�ƥ���˲�
%�汾ʱ�䣺2021/1/28 19:43
clc
clear
close all
load MandS_R2_u5.mat
gvar;    % ����ȫ�ֱ���
g2d=180/pi;
d2g=1/g2d;
att0=MandS.TR.att0;vn0=MandS.TR.vn0;pos0=MandS.TR.pos0;
ts=MandS.MINS.ts;
SINS=MandS.SINS;MINS=MandS.MINS;
TR=MandS.TR;
nn=2;nts=nn*ts;
len=length(SINS{1,2}.wis);
R=[2;0;0];
%% ���ߴ��ߵ�����
MERR.eb=[0;0;0]*dph;
MERR.web=[0;0;0]*dpsh;
MERR.db=[0;0;0]*ug;
MERR.wdb=[0;0;0]*ugpsHz;
%% �ӹ߸��崿�ߵ�����
% SERR.eb=[0;0;0]*dph;
% SERR.web=[0;0;0]*dpsh;
% SERR.db=[0;0;0]*ug;
% SERR.wdb=[0;0;0]*ugpsHz;
SERR.eb=[1;1;1]*dph;
SERR.web=[1;1;1]*dpsh;
SERR.db=[200;200;200]*ug;
SERR.wdb=[200;200;200]*ugpsHz;
%% �ٶ�ƥ�� ������
SINS_Ret_VL=zeros(floor(len/2),22);kk=1;t=0;
qnb_vm=a2qua(att0);
vn_vm=vn0;
posMINS_vm=pos0;
qnb_vs=a2qua(att0);
vn_vs=vn0;
posSINS_vs=pos0;
cl = cos(posSINS_vs(1,1)); Re = 6378137;
posSINS_vs(2,1)=posSINS_vs(2,1)+R(1)/(cl*Re);
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [ws1,vs1]=imuadderr(SINS{1,1}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,1}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL]=...
        my_Vrelinsupdate...
        (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
    SINS_Ret_VL(kk,:)=...
        [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
    kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
%% ��ͼ
msplot(211, SINS_Ret_VL(:,22), (SINS_Ret_VL(:,1:3)-SINS_Ret_VL(:,10:12))/arcdeg, 'Att / ( �� )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_VL(:,22), SINS_Ret_VL(:,4:6)+SINS_Ret_VL(:,19:21)-SINS_Ret_VL(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(211, SINS_Ret_VL(:,22), (SINS_Ret_VL(:,1:3)-SINS_Ret_VL(:,10:12))/arcdeg, 'Att / ( �� )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_VL(:,22), SINS_Ret_VL(:,4:6)-SINS_Ret_VL(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%% �ٶ�ƥ�� ������
SINS_Ret_Vu=zeros(floor(len/2),22);kk=1;t=0;
qnb_vm=a2qua(att0);
vn_vm=vn0;
posMINS_vm=pos0;
qnb_vs=a2qua(att0);
vn_vs=vn0;
posSINS_vs=pos0;
cl = cos(posSINS_vs(1,1)); Re = 6378137;
posSINS_vs(2,1)=posSINS_vs(2,1)+R(1)/(cl*Re);
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [ws1,vs1]=imuadderr(SINS{1,2}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,2}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL]=...
        my_Vrelinsupdate...
        (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
    SINS_Ret_Vu(kk,:)=...
        [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
    kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
%% ��ͼ
msplot(211, SINS_Ret_Vu(:,22), (SINS_Ret_Vu(:,1:3)-SINS_Ret_Vu(:,10:12))/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_Vu(:,22), SINS_Ret_Vu(:,4:6)+SINS_Ret_Vu(:,19:21)-SINS_Ret_Vu(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%% �ٶ�ƥ�� ������
SERR.eb=[1;1;1]*dph;
SERR.web=[1;1;1]*dpsh;
SERR.db=[200;200;200]*ug;
SERR.wdb=[200;200;200]*ugpsHz;
atterr0=[0;0;0]*arcdeg;
SINS_Ret_VLF=zeros(floor(len/2),22);kk=1;t=0;
qnb_vm=a2qua(att0);
vn_vm=vn0;
posMINS_vm=pos0;
qnb_vs=a2qua(att0+[1;2;3]*arcdeg);
vn_vs=vn0;
posSINS_vs=pos0;
cl = cos(posSINS_vs(1,1)); Re = 6378137;
posSINS_vs(2,1)=posSINS_vs(2,1)+R(1)/(cl*Re);
Filter.X=zeros(12,len/2);
Cmserr_1=my_a2mat(att0+atterr0);%ע���м䴦����ת˳��yxz ��׼ȷ�ĳ�ʼ��
qns=m2qua(Cmserr_1);%��ʼ��̬ ��׼ȷ�ĳ�ʼ��
fs0=[0;0;0];
KFinit.Qk = diag([SERR.wdb; SERR.web;])^2*nts;
KFinit.rk = [0.001;0.001;0.001];  
KFinit.Rk = diag(KFinit.rk)^2;
KFinit.P0 = diag([[0.001;0.001;0.001];[0.1;0.1;0.1]*arcdeg;SERR.eb; SERR.db])^2;
eth0 = earth(posMINS_vm, vn_vm);  % ������ز�������
kfft=my_kfftV12(eth0,q2mat(qns)',fs0,nts);
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
    [ws1,vs1]=imuadderr(SINS{1,1}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,1}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_vm,vn_vm,posMINS_vm,qnb_vs,vn_vs,posSINS_vs,VL,ethm]=...
        my_Vrelinsupdate...
        (qnb_vm,vn_vm,posMINS_vm,wm1,vm1,qnb_vs,vn_vs,posSINS_vs,ws1,vs1,R,ts);
    
    fs=sum(SINS{1,1}.fs(:,(2*k-1):(2*k)),2)/2;
    kfft=my_kfftV12(ethm,q2mat(qnb_vs)',fs,nts);
    kf.Phikk_1=kfft.phi;
    kf.Gammak=kfft.Gammak;
    
    ZV=vn_vs-(vn_vm+VL);
    
    kf = kfupdate(kf,ZV,'B');%�������˲�
    
    Filter.X(:,k)=kf.Xk;
%     qnb_vs = qdelphi(qnb_vs,-kf.Xk(4:6));  kf.Xk(4:6) = 0;  % ������
%     vn_vs = vn_vs - kf.Xk(1:3);  kf.Xk(1:3) = 0;
    
    SINS_Ret_VLF(kk,:)=...
        [q2att(qnb_vm);vn_vm;posMINS_vm;q2att(qnb_vs);vn_vs;posSINS_vs;VL;t];
    kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
Filter.XT=Filter.X';
%% ��ͼ
msplot(211, SINS_Ret_VLF(:,22), (SINS_Ret_VLF(:,1:3)-SINS_Ret_VLF(:,10:12))/arcdeg, 'Att / ( �� )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)+SINS_Ret_VLF(:,19:21)-SINS_Ret_VLF(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(211, SINS_Ret_VLF(:,22), (SINS_Ret_VLF(:,1:3)-SINS_Ret_VLF(:,10:12))/arcdeg, 'Att / ( �� )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(212, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)-SINS_Ret_VLF(:,13:15), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%�Ա��˲������ٶ����
msplot(311, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)+SINS_Ret_VLF(:,19:21), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')%���������ߵ����
msplot(312, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,13:15)-Filter.XT(:,1:3), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(313, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,4:6)+SINS_Ret_VLF(:,19:21)-SINS_Ret_VLF(:,13:15)+Filter.XT(:,1:3), 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
%% �Ա��˲�������̬���
msplot(311, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,1:3)/arcdeg,'Att / ( �� )'); legend('\it\theta', '\it\gamma', '\it\psi')%���������ߵ����
msplot(312, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,10:12)/arcdeg,'Att / (�� )'); legend('\it\theta', '\it\gamma', '\it\psi')%���������ߵ����
msplot(313, SINS_Ret_VLF(:,22), SINS_Ret_VLF(:,1:3)/arcdeg-SINS_Ret_VLF(:,10:12)/arcdeg,'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')