function [INSSFR,Filter]=my_getFResult(MINS,SINS,KFinit,atterr0,looplen,varargin)
% �������ƣ�my_SINSgetFResult
% �������ܣ����ӹ����������Ե����˲�
% ���룺U      :��ǰ����ϵ������ٶ�
%       varargin:1�������ʾSERR
%                2�������ʾws_m_addnoise fs_m_addnoise
%                                
%                                
%                               
% �����INSSFR:��Ե������
gvar;    % ����ȫ�ֱ���
looplen=floor(looplen);
Rloop=SINS.R0;
Uloop=cross(MINS.wim0,SINS.R0);
ts=SINS.ts;
nts=2*ts;
INSSFR.Rall=zeros(looplen,3);
INSSFR.vnall=zeros(looplen,3);
INSSFR.qmsall=zeros(looplen,4);
INSSFR.attall=zeros(3,looplen);

Filter.X=zeros(15,looplen);

Cmserr_1=my_a2mat(atterr0);%ע���м䴦����ת˳��yxz ��׼ȷ�ĳ�ʼ��
qms=m2qua(Cmserr_1);%��ʼ��̬ ��׼ȷ�ĳ�ʼ��
attint=m2att(Cmserr_1);%��ת����yxz ������θ����� ����� ��λ��

fs0=[0;0;0];
kfft=my_kfft15(MINS.wim0,q2mat(qms)',fs0,nts);
kf = kfinit(KFinit.Qk, KFinit.Rk, KFinit.P0,kfft.phi,kfft.H);  % kf�˲�����ʼ��
kf.Xk(10:12) = [0.1;0.1;0.1]*dph;
kf.Xk(13:15) = [20;20;20]*ug;
RFilter=SINS.R0;
for i=1:looplen
    if nargin==6
        SERR=varargin{1};
        [ws_m, fs_m] = my_imuadderr(SINS.wis(:,(2*i-1):(2*i))',...
            SINS.fs(:,(2*i-1):(2*i))', ...
            SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);%�ӹ�ע������
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; fm_m=MINS.fm(:,(2*i-1):(2*i))';
    elseif nargin==7
        ws_m=varargin{1}(:,:,i);
        fs_m=varargin{2}(:,:,i);
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; fm_m=MINS.fm(:,(2*i-1):(2*i))';
    else
        ws_m=SINS.wis(:,(2*i-1):(2*i))';
        fs_m=SINS.fs(:,(2*i-1):(2*i))';
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; fm_m=MINS.fm(:,(2*i-1):(2*i))';
    end


    [qms,vn,Uloop,Rloop] = my_relinsupdate5(qms,Uloop,Rloop,wm_m, fm_m, ws_m, fs_m, ts);% �ӹߴ���Ե���
    
    RFilter=my_Ofilter(RFilter,Rloop,2,nts);%��ͨ�˲�
    
    attnow=q2att(qms);
    uf=attnow-attint;
    
    dR=my_getdR(RFilter,SINS.R0,uf);%�������
    
%     wim=1.5*wm_m(2,:)-0.5*wm_m(1,:);%����������
%     fs=1.5*fm_m(2,:)-0.5*fm_m(1,:);
    wim=sum(wm_m)/2;%����������
    fs=sum(fm_m)/2;
    kfft=my_kfft15(wim',q2mat(qms)',fs',nts);
    kf.Phikk_1=kfft.phi;
    kf.Gammak=kfft.Gammak;
    
    kf = kfupdate(kf,dR,'B');%�������˲�
    
    Filter.X(:,i)=kf.Xk;
    attint=q2att(qdelphi(a2qua(attint),-kf.Xk(1:3)));% 2021 03 15���� ���³�ʼ��̬�����
%     [rz,rx,ry]=dcm2angle(q2mat(rv2q(kf.Xk(1:3))),'zxy');
%     Filter.atterr=[rx;ry;-rz];
    
%     kf.Xk(1:3) = 0;
    qms = qdelphi(qms,-kf.Xk(1:3));  kf.Xk(1:3) = 0;  % ������
    Uloop = Uloop - kf.Xk(4:6);  kf.Xk(4:6) = 0;
    Rloop = Rloop - kf.Xk(7:9);  kf.Xk(7:9) = 0;
    
    %����
    qint = qdelphi(m2qua(a2mat(attint)),-kf.Xk(1:3));%������ת��ʼ��̬
    attint=q2att(qint);
    
    INSSFR.qmsall(i,:)=qms';
    INSSFR.Rall(i,:)=Rloop';
    INSSFR.vnall(i,:)=vn';
    
    [rz,rx,ry]=dcm2angle(q2mat(INSSFR.qmsall(i,:)),'zxy');%ת������ͳ������
    INSSFR.attall(:,i)=[rx,ry,-rz];
    
    if mod(i,looplen/100)==0
              fprintf('\b\b\b\b\b\b\b%5.0f %%', i/looplen*100);			%������ʾ
    end
end

% for i=1:length(INSSFR.qmsall)%���ӹ߽��������̬�Ƕ�
%     [rz,rx,ry]=dcm2angle(q2mat(INSSFR.qmsall(i,:)),'zxy');%ת������ͳ������
%     INSSFR.attall(:,i)=[rx,ry,-rz];
% end