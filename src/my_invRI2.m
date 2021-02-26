function INSS=my_invRI2(R,att,attm,atterr0,wim,fm,ts,varargin)
% 函数名称：my_invRI
% 函数功能：反解算相对导航子惯无噪声数据
% 输入：U      :当前导航系下相对速度
%       varargin:函数除微分项其余参数,以元组给入          
%              每一个参数中依次存放Cms，（主到子坐标变换阵，与Csm区分）3*3
%                                  dU  速度微分
%                                  fm  主惯比力
%                                  wim 主惯相对于惯性坐标系速度在主惯坐标系的投影
% 输出：dU:当前导航系下相对速度微分
R0=varargin{1};
att0=varargin{2};
U0=varargin{3};%U0=diffR0+cross(wim0,R0);
len=size(R,2);

diffR=my_diff(R,2,R0,1)/ts;

attn=zeros(3,len);
for i=1:len
    attn(:,i)=m2att(a2mat(-att(:,i))*a2mat(attm(:,i)));
end

wis_m=att2wm(attn,m2att(a2mat(-att0)*a2mat(attm(:,1))))/ts;%主惯坐标系下主子相对角速度
% wis_m=wim+wms_m;%主惯坐标系下子惯相对于惯性系角速度
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
diffU=my_diff(INSS.U,2,U0,1)/ts;
for i=1:len
    INSS.fs(:,i)=my_dUtofs(INSS.U(:,i),Cms_cell{i},diffU(:,i),fm(:,i),wim(:,i));
end
INSS.Cms=Cms_cell;
