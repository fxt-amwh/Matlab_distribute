function INSS=my_invRI(R,att,atterr0,wim,fm,ts,varargin)
% �������ƣ�my_invRI
% �������ܣ���������Ե����ӹ�����������
% ���룺1 R                 :�ӹ��������е��˶�λ�þ���
%       2 att               :������΢�����������,��Ԫ�����          
%       3 atterr0           :�ӹ߳�ʼ��װ���
%       4 wim               :���߽��ٶȾ���
%       5 fm                :���߱�������
%       6 ts                :����ʱ��
%       7 varargin{1}       :R0
%       8 varargin{2}       :att0
%       9 varargin{3}       :U0
% �����INSS:�ӹ߽ṹ��
%% ������ʼ��
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

wms_m=-att2wm(att,att0)/ts;%��������ϵ��������Խ��ٶ�
wis_m=wim+wms_m;%��������ϵ���ӹ�����ڹ���ϵ���ٶ�
INSS.wis=zeros(3,len);

Cms_cell=cell(1,len);
Cmserr=my_a2mat(atterr0);%ע���м䴦����ת˳��yxz

INSS.fs=zeros(3,len);
INSS.U=zeros(3,len);

for i=1:len
    Cms_cell{i}=Cmserr*a2mat(att(:,i));%��������任����
%     Cms_cell{i}=a2mat(att(:,i))*Cmserr;%��������任����
    INSS.wis(:,i)=Cms_cell{i}*wis_m(:,i);
    INSS.U(:,i)=my_dRtoU(R(:,i),diffR(:,i),wim(:,i));
    
end
diffU=my_diff(INSS.U,2,U0,2)/ts;
for i=1:len
    INSS.fs(:,i)=my_dUtofs(INSS.U(:,i),Cms_cell{i},diffU(:,i),fm(:,i),wim(:,i));
end
%% �������
INSS.Cms=Cms_cell;
INSS.atttrue=zeros(3,len);
INSS.attuftrue=zeros(3,len);
for i=1:len%��ʵ�Ƕ�
   [rz,rx,ry]=dcm2angle(INSS.Cms{i},'zxy');
   INSS.atttrue(:,i)=[rx,ry,-rz];   %�а�װ���ǵ���̬��ֵ
   [rz,rx,ry]=dcm2angle(a2mat(att(:,i)),'zxy');
   INSS.attuftrue(:,i)=[rx,ry,-rz]; %û�а�װ���ǵ���̬��ֵ
end
INSS.R=R;
INSS.R0=R0;
INSS.Cmserr=Cmserr;
INSS.ts=ts;