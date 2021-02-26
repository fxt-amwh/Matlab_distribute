%�ļ����ƣ�test1
%������Ϣ��������
%�������������ӹߵ��ֱ���е����õ��ٶ�
%�汾ʱ�䣺2021/1/27 17:19
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
%% ���ߴ��ߵ�����
MINS_Ret=zeros(floor(len/2),10);kk=1;t=0;
MERR.eb=[0;0;0]*dph;
MERR.web=[0;0;0]*dpsh;
MERR.db=[0;0;0]*ug;
MERR.wdb=[0;0;0]*ugpsHz;
qnb=a2qua(att0);
vn=vn0;
posMINS=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(MINS.wim(:,(2*k-1):(2*k))'*ts, ...
        MINS.fm(:,(2*k-1):(2*k))'*ts, ...
        MERR.eb, MERR.web, MERR.db, MERR.wdb, ts);
    [qnb,vn,posMINS]=my_insupdate(qnb,vn,posMINS,wm1,vm1,ts);
    MINS_Ret(kk,:)=[q2att(qnb);vn;posMINS;t];kk=kk+1;
    
    if mod(t,100)<nts,disp(fix(t));end
end
%% �ӹ߸��崿�ߵ�����
SINS_Ret_L=zeros(floor(len/2),10);kk=1;t=0;
SERR.eb=[0;0;0]*dph;
SERR.web=[0;0;0]*dpsh;
SERR.db=[0;0;0]*ug;
SERR.wdb=[0;0;0]*ugpsHz;
qnb_L=a2qua(att0);
vn_L=vn0;
posSINS_L=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(SINS{1,1}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,1}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_L,vn_L,posSINS_L]=my_insupdate(qnb_L,vn_L,posSINS_L,wm1,vm1,ts);
    SINS_Ret_L(kk,:)=[q2att(qnb_L);vn_L;posSINS_L;t];kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
%% �ӹ��������ߵ�����
SINS_Ret_u=zeros(floor(len/2),10);kk=1;t=0;
qnb_u=a2qua(att0);
vn_u=vn0;
posSINS_u=pos0;
for k=1:len/2
    t=t+nts;
    [wm1,vm1]=imuadderr(SINS{1,2}.wis(:,(2*k-1):(2*k))'*ts, ...
        SINS{1,2}.fs(:,(2*k-1):(2*k))'*ts, ...
        SERR.eb, SERR.web, SERR.db, SERR.wdb, ts);
    [qnb_u,vn_u,posSINS_u]=my_insupdate(qnb_u,vn_u,posSINS_u,wm1,vm1,ts);
    SINS_Ret_u(kk,:)=[q2att(qnb_u);vn_L;posSINS_u;t];kk=kk+1;
    if mod(t,100)<nts,disp(fix(t));end
end
%% ��ͼ
delpos=deltapos(TR.pos);
delposMINS=deltapos(MINS_Ret(:,7:9));
delposSINS_L=deltapos(SINS_Ret_L(:,7:9));
delposSINS_u=deltapos(SINS_Ret_u(:,7:9));
%��ʵ�켣�����߶Ա�
figure 
plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
hold on;
plot3(delposMINS(:,1),delposMINS(:,2),delposMINS(:,3));
grid on;
xlabel('γ�Ⱦ���/m'),ylabel('���Ⱦ���/m'),zlabel('��/m');
err=delposMINS(end,:)-delpos(end,:);
disp 'γ����� ������� �߶����'
disp(err)
%��ʵ�켣���ӹ߸���Ա�
figure
plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
hold on;
plot3(delposSINS_L(:,1),delposSINS_L(:,2),delposSINS_L(:,3));
grid on;
xlabel('γ�Ⱦ���/m'),ylabel('���Ⱦ���/m'),zlabel('��/m');
err=delposSINS_L(end,:)-delpos(end,:);
disp 'γ����� ������� �߶����'
disp(err)
%��ʵ�켣���ӹ������Ա�
figure
plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
hold on;
plot3(delposSINS_u(:,1),delposSINS_u(:,2),delposSINS_u(:,3));
grid on;
xlabel('γ�Ⱦ���/m'),ylabel('���Ⱦ���/m'),zlabel('��/m');
err=delposSINS_u(end,:)-delpos(end,:);
disp 'γ����� ������� �߶����'
disp(err)
%�ӹ߸������ӹ������Ա�
figure
plot3(delposSINS_L(:,1),delposSINS_L(:,2),delposSINS_L(:,3),'LineWidth',5);
hold on;
plot3(delposSINS_u(:,1),delposSINS_u(:,2),delposSINS_u(:,3));
grid on;
xlabel('γ�Ⱦ���/m'),ylabel('���Ⱦ���/m'),zlabel('��/m');
err=delposSINS_u(end,:)-delposSINS_L(end,:);
disp 'γ����� ������� �߶����'
disp(err)
%��ʵ�켣���������
figure 
subplot(311)
plot(MINS_Ret(:,10),delpos(1:2:end,1)-delposMINS(:,1));ylabel('γ�Ⱦ���/m')
title("��ʵ�켣������γ�Ⱦ������");grid on;
subplot(312)
plot(MINS_Ret(:,10),delpos(1:2:end,2)-delposMINS(:,2));ylabel('���Ⱦ���/m')
title("��ʵ�켣�����߾��Ⱦ������");grid on;
subplot(313)
plot(MINS_Ret(:,10),delpos(1:2:end,3)-delposMINS(:,3));ylabel('��/m');
grid on;
title("��ʵ�켣�����߸߶����");
% %��ʵ�켣���ӹ߸������
% figure 
% subplot(311)
% plot(MINS_Ret(:,10),delpos(1:2:end,1)-delposSINS_L(:,1));ylabel('γ�Ⱦ���/m')
% title("��ʵ�켣���ӹ߸���γ�Ⱦ������");grid on;
% subplot(312)
% plot(MINS_Ret(:,10),delpos(1:2:end,2)-delposSINS_L(:,2));ylabel('���Ⱦ���/m')
% title("��ʵ�켣���ӹ߸��徭�Ⱦ������");grid on;
% subplot(313)
% plot(MINS_Ret(:,10),delpos(1:2:end,3)-delposSINS_L(:,3));ylabel('��/m');
% grid on;
% title("��ʵ�켣���ӹ߸���߶����");
% %��ʵ�켣���ӹ��������
% figure 
% subplot(311)
% plot(MINS_Ret(:,10),delpos(1:2:end,1)-delposSINS_u(:,1));ylabel('γ�Ⱦ���/m')
% title("��ʵ�켣���ӹ�����γ�Ⱦ������");grid on;
% subplot(312)
% plot(MINS_Ret(:,10),delpos(1:2:end,2)-delposSINS_u(:,2));ylabel('���Ⱦ���/m')
% title("��ʵ�켣���ӹ��������Ⱦ������");grid on;
% subplot(313)
% plot(MINS_Ret(:,10),delpos(1:2:end,3)-delposSINS_u(:,3));ylabel('��/m');
% grid on;
% title("��ʵ�켣���ӹ������߶����");
%���߹켣���ӹ߸������
figure 
subplot(311)
plot(MINS_Ret(:,10),delposMINS(:,1)-delposSINS_L(:,1));ylabel('γ�Ⱦ���/m')
title("���߹켣���ӹ߸���γ�Ⱦ������");grid on;
subplot(312)
plot(MINS_Ret(:,10),delposMINS(:,2)-delposSINS_L(:,2));ylabel('���Ⱦ���/m')
title("���߹켣���ӹ߸��徭�Ⱦ������");grid on;
subplot(313)
plot(MINS_Ret(:,10),delposMINS(:,3)-delposSINS_L(:,3));ylabel('��/m');
grid on;
title("���߹켣���ӹ߸���߶����");
%��ʵ�켣���ӹ��������
figure 
subplot(311)
plot(MINS_Ret(:,10),delposMINS(:,1)-delposSINS_u(:,1));ylabel('γ�Ⱦ���/m')
title("���߹켣���ӹ�����γ�Ⱦ������");grid on;
subplot(312)
plot(MINS_Ret(:,10),delposMINS(:,2)-delposSINS_u(:,2));ylabel('���Ⱦ���/m')
title("���߹켣���ӹ��������Ⱦ������");grid on;
subplot(313)
plot(MINS_Ret(:,10),delposMINS(:,3)-delposSINS_u(:,3));ylabel('��/m');
grid on;
title("���߹켣���ӹ������߶����");
%�ӹ߸������ӹ��������
figure 
subplot(311)
plot(MINS_Ret(:,10),delposSINS_L(:,1)-delposSINS_u(:,1));ylabel('γ�Ⱦ���/m')
title("�ӹ߸������ӹ�����γ�Ⱦ������");grid on;
subplot(312)
plot(MINS_Ret(:,10),delposSINS_L(:,2)-delposSINS_u(:,2));ylabel('���Ⱦ���/m')
title("�ӹ߸������ӹ��������Ⱦ������");grid on;
subplot(313)
plot(MINS_Ret(:,10),delposSINS_L(:,3)-delposSINS_u(:,3));ylabel('��/m');
grid on;
title("�ӹ߸������ӹ������߶����");