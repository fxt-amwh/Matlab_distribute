function INSSR=my_SINSgetResult(MINS,SINS,looplen,varargin)
% 函数名称：my_SINSgetResult
% 函数功能：由子惯输出进行相对导航
% 输入：MINS     ：主惯结构体
%       SINS     ：子惯结构体
%       varargin:函数除微分项其余参数,以元组给入          
%              每一个参数中依次存放Cms，（主到子坐标变换阵，与Csm区分）3*3
%                                  dU  速度微分
%                                  fm  主惯比力
%                                  wim 主惯相对于惯性坐标系速度在主惯坐标系的投影
% 输出：INSSR:相对导航结果
%% 相对导航
looplen=floor(looplen);
SERR=varargin{1};
ts=SINS.ts;
Rloop=SINS.R0;
Uloop=cross(MINS.wim0,SINS.R0);
if nargin==4
    qmsloop=m2qua(SINS.Cmserr);%初始相对姿态对应的四元数
else
    atterr0=varargin{2};
    Cmserr=my_a2mat(atterr0);%注意中间处理旋转顺序yxz
    qmsloop=m2qua(Cmserr);%初始相对姿态对应的四元数
end



%相对导航结果变量初始化
INSSR.Rall=zeros(looplen,3);%相对位置
INSSR.vnall=zeros(looplen,3);%相对速度
INSSR.qmsall=zeros(looplen,4);%相对姿态四元数
INSSR.attall=zeros(3,looplen);%相对姿态 旋转顺序zxy
%子惯注入噪声结果变量初始化
INSSR.ws_m_addnoise=zeros(2,3,looplen);
INSSR.fs_m_addnoise=zeros(2,3,looplen);
for i=1:looplen
    if mod(i,looplen/100)==0
              fprintf('\b\b\b\b\b\b\b%5.0f %%', i/looplen*100);			%进度显示
    end
    [ws_m, fs_m] = my_imuadderr(SINS.wis(:,(2*i-1):(2*i))', SINS.fs(:,(2*i-1):(2*i))', SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);%子惯注入噪声
    INSSR.ws_m_addnoise(:,:,i)=ws_m;%保存注入噪声的子惯数据
    INSSR.fs_m_addnoise(:,:,i)=fs_m;%保存注入噪声的子惯数据
    wm_m=MINS.wim(:,(2*i-1):(2*i))'; fm_m=MINS.fm(:,(2*i-1):(2*i))';%获取主惯数据
    [qmsloop,vn,Uloop,Rloop] = my_relinsupdate5(qmsloop,Uloop,Rloop,wm_m, fm_m, ws_m, fs_m, ts);% 子惯纯相对导航
    %保存相对解算结果
    INSSR.qmsall(i,:)=qmsloop';
    INSSR.Rall(i,:)=Rloop';
    INSSR.vnall(i,:)=vn';
    [rz,rx,ry]=dcm2angle(q2mat(INSSR.qmsall(i,:)),'zxy');%转换到传统的坐标
    INSSR.attall(:,i)=[rx,ry,-rz];
end
% for i=1:length(INSSR.qmsall)%求子惯解算出的姿态角度
%     [rz,rx,ry]=dcm2angle(q2mat(INSSR.qmsall(i,:)),'zxy');%转换到传统的坐标
%     INSSR.attall(:,i)=[rx,ry,-rz];
% end