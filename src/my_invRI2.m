function INSS=my_invRI2(R,att,attm,atterr0,wim,fm,ts,varargin)
% �������ƣ�my_invRI
% �������ܣ���������Ե����ӹ�����������
% ���룺U      :��ǰ����ϵ������ٶ�
%       varargin:������΢�����������,��Ԫ�����          
%              ÿһ�����������δ��Cms��������������任����Csm���֣�3*3
%                                  dU  �ٶ�΢��
%                                  fm  ���߱���
%                                  wim ��������ڹ�������ϵ�ٶ�����������ϵ��ͶӰ
% �����dU:��ǰ����ϵ������ٶ�΢��
R0=varargin{1};
att0=varargin{2};
U0=varargin{3};%U0=diffR0+cross(wim0,R0);
len=size(R,2);

diffR=my_diff(R,2,R0,1)/ts;

attn=zeros(3,len);
for i=1:len
    attn(:,i)=m2att(a2mat(-att(:,i))*a2mat(attm(:,i)));
end

wis_m=att2wm(attn,m2att(a2mat(-att0)*a2mat(attm(:,1))))/ts;%��������ϵ��������Խ��ٶ�
% wis_m=wim+wms_m;%��������ϵ���ӹ�����ڹ���ϵ���ٶ�
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
diffU=my_diff(INSS.U,2,U0,1)/ts;
for i=1:len
    INSS.fs(:,i)=my_dUtofs(INSS.U(:,i),Cms_cell{i},diffU(:,i),fm(:,i),wim(:,i));
end
INSS.Cms=Cms_cell;
