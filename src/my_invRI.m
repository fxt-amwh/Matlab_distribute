function INSS=my_invRI(R,att,atterr0,wim,fm,ts,varargin)
% 函数名称：my_invRI
% 函数功能：反解算相对导航子惯无噪声数据
% 输入：1 R                 :子惯在主惯中的运动位置矩阵
%       2 att               :函数除微分项其余参数,以元组给入          
%       3 atterr0           :子惯初始安装误差
%       4 wim               :主惯角速度矩阵
%       5 fm                :主惯比力矩阵
%       6 ts                :采样时间
%       7 varargin{1}       :R0
%       8 varargin{2}       :att0
%       9 varargin{3}       :U0
% 输出：INSS:子惯结构体
%% 变量初始化
R0=R(:,1);
att0=att(:,1);
U0=[0;0;0]+cross(wim(:,1),R0);
if nargin==7
    R0=varargin{1};
elseif nargin==8
    R0=varargin{1};
    att0=varargin{2};
elseif nargin==9
    R0=varargin{1};
    att0=varargin{2};
    U0=varargin{3};%U0=diffR0+cross(wim0,R0);
end

len=size(R,2);

diffR=my_diff(R,2,R0,1)/ts;

wms_m=-att2wm(att,att0)/ts;%主惯坐标系下主子相对角速度
wis_m=wim+wms_m;%主惯坐标系下子惯相对于惯性系角速度
INSS.wis=zeros(3,len);

Cms_cell=cell(1,len);
Cmserr=my_a2mat(atterr0);%注意中间处理旋转顺序yxz

INSS.fs=zeros(3,len);
INSS.U=zeros(3,len);

for i=1:len
    Cms_cell{i}=Cmserr*a2mat(att(:,i));%主子坐标变换矩阵
%     Cms_cell{i}=a2mat(att(:,i))*Cmserr;%主子坐标变换矩阵
    INSS.wis(:,i)=Cms_cell{i}*wis_m(:,i);
    INSS.U(:,i)=my_dRtoU(R(:,i),diffR(:,i),wim(:,i));
    
end
diffU=my_diff(INSS.U,2,U0,2)/ts;
for i=1:len
    INSS.fs(:,i)=my_dUtofs(INSS.U(:,i),Cms_cell{i},diffU(:,i),fm(:,i),wim(:,i));
end
%% 输出配置
INSS.Cms=Cms_cell;
INSS.atttrue=zeros(3,len);
INSS.attuftrue=zeros(3,len);
for i=1:len%真实角度
   [rz,rx,ry]=dcm2angle(INSS.Cms{i},'zxy');
   INSS.atttrue(:,i)=[rx,ry,-rz];   %有安装误差角的姿态真值
   [rz,rx,ry]=dcm2angle(a2mat(att(:,i)),'zxy');
   INSS.attuftrue(:,i)=[rx,ry,-rz]; %没有安装误差角的姿态真值
end
INSS.R=R;
INSS.R0=R0;
INSS.Cmserr=Cmserr;
INSS.ts=ts;