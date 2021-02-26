function INSSR=my_SINSgetResult(MINS,SINS,looplen,varargin)
% �������ƣ�my_SINSgetResult
% �������ܣ����ӹ����������Ե���
% ���룺MINS     �����߽ṹ��
%       SINS     ���ӹ߽ṹ��
%       varargin:������΢�����������,��Ԫ�����          
%              ÿһ�����������δ��Cms��������������任����Csm���֣�3*3
%                                  dU  �ٶ�΢��
%                                  fm  ���߱���
%                                  wim ��������ڹ�������ϵ�ٶ�����������ϵ��ͶӰ
% �����INSSR:��Ե������
%% ��Ե���
looplen=floor(looplen);
SERR=varargin{1};
ts=SINS.ts;
Rloop=SINS.R0;
Uloop=cross(MINS.wim0,SINS.R0);
if nargin==4
    qmsloop=m2qua(SINS.Cmserr);%��ʼ�����̬��Ӧ����Ԫ��
else
    atterr0=varargin{2};
    Cmserr=my_a2mat(atterr0);%ע���м䴦����ת˳��yxz
    qmsloop=m2qua(Cmserr);%��ʼ�����̬��Ӧ����Ԫ��
end



%��Ե������������ʼ��
INSSR.Rall=zeros(looplen,3);%���λ��
INSSR.vnall=zeros(looplen,3);%����ٶ�
INSSR.qmsall=zeros(looplen,4);%�����̬��Ԫ��
INSSR.attall=zeros(3,looplen);%�����̬ ��ת˳��zxy
%�ӹ�ע���������������ʼ��
INSSR.ws_m_addnoise=zeros(2,3,looplen);
INSSR.fs_m_addnoise=zeros(2,3,looplen);
for i=1:looplen
    if mod(i,looplen/100)==0
              fprintf('\b\b\b\b\b\b\b%5.0f %%', i/looplen*100);			%������ʾ
    end
    [ws_m, fs_m] = my_imuadderr(SINS.wis(:,(2*i-1):(2*i))', SINS.fs(:,(2*i-1):(2*i))', SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);%�ӹ�ע������
    INSSR.ws_m_addnoise(:,:,i)=ws_m;%����ע���������ӹ�����
    INSSR.fs_m_addnoise(:,:,i)=fs_m;%����ע���������ӹ�����
    wm_m=MINS.wim(:,(2*i-1):(2*i))'; fm_m=MINS.fm(:,(2*i-1):(2*i))';%��ȡ��������
    [qmsloop,vn,Uloop,Rloop] = my_relinsupdate5(qmsloop,Uloop,Rloop,wm_m, fm_m, ws_m, fs_m, ts);% �ӹߴ���Ե���
    %������Խ�����
    INSSR.qmsall(i,:)=qmsloop';
    INSSR.Rall(i,:)=Rloop';
    INSSR.vnall(i,:)=vn';
    [rz,rx,ry]=dcm2angle(q2mat(INSSR.qmsall(i,:)),'zxy');%ת������ͳ������
    INSSR.attall(:,i)=[rx,ry,-rz];
end
% for i=1:length(INSSR.qmsall)%���ӹ߽��������̬�Ƕ�
%     [rz,rx,ry]=dcm2angle(q2mat(INSSR.qmsall(i,:)),'zxy');%ת������ͳ������
%     INSSR.attall(:,i)=[rx,ry,-rz];
% end