function [INSSFR,Filter]=my_getFResult(MINS,SINS,KFinit,atterr0,looplen,varargin)
% 函数名称：my_SINSgetFResult
% 函数功能：由子惯输出进行相对导航滤波
% 输入：U      :当前导航系下相对速度
%       varargin:1个输入表示SERR
%                2个输入表示ws_m_addnoise fs_m_addnoise
%                                
%                                
%                               
% 输出：INSSFR:相对导航结果
gvar;    % 加载全局变量
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

Cmserr_1=my_a2mat(atterr0);%注意中间处理旋转顺序yxz 不准确的初始角
qms=m2qua(Cmserr_1);%初始姿态 不准确的初始角
attint=m2att(Cmserr_1);%旋转次序yxz 输出依次俯仰角 横滚角 方位角

fs0=[0;0;0];
kfft=my_kfft15(MINS.wim0,q2mat(qms)',fs0,nts);
kf = kfinit(KFinit.Qk, KFinit.Rk, KFinit.P0,kfft.phi,kfft.H);  % kf滤波器初始化
kf.Xk(10:12) = [0.1;0.1;0.1]*dph;
kf.Xk(13:15) = [20;20;20]*ug;
RFilter=SINS.R0;
for i=1:looplen
    if nargin==6
        SERR=varargin{1};
        [ws_m, fs_m] = my_imuadderr(SINS.wis(:,(2*i-1):(2*i))',...
            SINS.fs(:,(2*i-1):(2*i))', ...
            SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);%子惯注入噪声
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


    [qms,vn,Uloop,Rloop] = my_relinsupdate5(qms,Uloop,Rloop,wm_m, fm_m, ws_m, fs_m, ts);% 子惯纯相对导航
    
    RFilter=my_Ofilter(RFilter,Rloop,2,nts);%低通滤波
    
    attnow=q2att(qms);
    uf=attnow-attint;
    
    dR=my_getdR(RFilter,SINS.R0,uf);%求出量测
    
%     wim=1.5*wm_m(2,:)-0.5*wm_m(1,:);%可能有问题
%     fs=1.5*fm_m(2,:)-0.5*fm_m(1,:);
    wim=sum(wm_m)/2;%可能有问题
    fs=sum(fm_m)/2;
    kfft=my_kfft15(wim',q2mat(qms)',fs',nts);
    kf.Phikk_1=kfft.phi;
    kf.Gammak=kfft.Gammak;
    
    kf = kfupdate(kf,dR,'B');%卡尔曼滤波
    
    Filter.X(:,i)=kf.Xk;
    attint=q2att(qdelphi(a2qua(attint),-kf.Xk(1:3)));% 2021 03 15新增 更新初始姿态误差阵
%     [rz,rx,ry]=dcm2angle(q2mat(rv2q(kf.Xk(1:3))),'zxy');
%     Filter.atterr=[rx;ry;-rz];
    
%     kf.Xk(1:3) = 0;
    qms = qdelphi(qms,-kf.Xk(1:3));  kf.Xk(1:3) = 0;  % ·反馈
    Uloop = Uloop - kf.Xk(4:6);  kf.Xk(4:6) = 0;
    Rloop = Rloop - kf.Xk(7:9);  kf.Xk(7:9) = 0;
    
    %新增
    qint = qdelphi(m2qua(a2mat(attint)),-kf.Xk(1:3));%尝试旋转初始姿态
    attint=q2att(qint);
    
    INSSFR.qmsall(i,:)=qms';
    INSSFR.Rall(i,:)=Rloop';
    INSSFR.vnall(i,:)=vn';
    
    [rz,rx,ry]=dcm2angle(q2mat(INSSFR.qmsall(i,:)),'zxy');%转换到传统的坐标
    INSSFR.attall(:,i)=[rx,ry,-rz];
    
    if mod(i,looplen/100)==0
              fprintf('\b\b\b\b\b\b\b%5.0f %%', i/looplen*100);			%进度显示
    end
end

% for i=1:length(INSSFR.qmsall)%求子惯解算出的姿态角度
%     [rz,rx,ry]=dcm2angle(q2mat(INSSFR.qmsall(i,:)),'zxy');%转换到传统的坐标
%     INSSFR.attall(:,i)=[rx,ry,-rz];
% end