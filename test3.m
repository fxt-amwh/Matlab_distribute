%% һ��+2�ӣ��м䱸�ݣ�
% 1��������ƹ켣�˶����������ߴ��������
% 2��������ƻ����˶��������ӹߴ���������
% 3���ӹ�����ע�����
% 4����������Ե���
% 5���˲�
% 6��������
%% �����˶�������Ʒ���
clc
clear
close all
gvar;    % ����ȫ�ֱ���
ts = 0.1;%����ʱ��
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;105];
%     ���������� ��������� ��λ������ ������ٶ� ����ʱ�� �����ι켣��ƣ�
% wat = [ 0,         -60,        120,       0,         30       %��ֹ30
%         0,         60,         -120,      0,         60       %����90
%         0,         0,          0,         1,         20       %����110
%         60,        0,          0,         0,         30        %����140
%         0,         0,          0,         0,         10       %����150
%         -60,       0,          0,         0,         30        %����180
%         0,         -60,        120,       0,         60       %��ֹ240
%         0,         60,         -120,      0,         60       %����300
%         0,         0,          0,         0,         30       %����330
%         0,         0,          60,        -1,        20       %����350
%         60,        0,          0,         0,         30        %����380
%         0,         0,          0,         0,         10       %����390
%         -60,       0,          0,         0,         30        %����420
%         0,         0,          0,         0,         60       %����480
%         -60,       0,          0,         0,         30        %����510
%         0,         0,          0,         0,         10       %����520
%         60,        0,          0,         0,         30        %����550
%         0,         0,          0,         0,         60       %����610
%         0,         0,          60,        0,         20        %����630
%         0,         0,          -60,       0,         20        %����650
%         -60,       0,          0,         0,         30        %����680
%         0,         0,          0,         0,         10       %����690
%         60,        0,          0,         0,         30        %����720
%         
%         ];    %��ֹ
wat = [  0,         0,          0,         0,         10       %��ֹ10
        60,         0,          0,         0,         60       %����70
        -60,        0,          0,         0,         60       %����130
        0,          0,          0,         0,         50       %����180
        0,          0,          60,       0,         360      %����ת��540
        0,          0,          0,         0,         50       %����590
        -60,        0,          0,         0,         60       %����650
        60,         0,          0,         0,        60       %����710
        0,          0,          0,         0,         10       %����720
        ];    %��ֹ
wat(:,1:3) = wat(:,1:3)*arcdeg/60;  % ->deg/s
fprintf('�켣������...'); 
[attm, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
fprintf('�켣������ɣ�\n���߷�����...');
[wm, vm] = my_av2imu(attm, vn, pos, ts);
tt = (0:length(attm)-1)'*ts;
% �켣��ͼ
msplot(221, tt, attm/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(222, tt, vn, 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(223, tt, deltapos(pos), '\DeltaPos / m');
      legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')
msplot(224, pos(:,2)/arcdeg, pos(:,1)/arcdeg, '\itL\rm / ( \circ )', '\it\lambda\rm / ( \circ)');
      hold on, plot(pos(1,2)/arcdeg, pos(1,1)/arcdeg, 'ro');
%  ���Դ�������Ϣ��ͼ
msplot(121, tt(2:end), wm/ts/arcdeg, '\it\omega^b_{ib}\rm / ( \circ.s^{-1} )');
      legend('\it\omega^b_{ibx}', '\it\omega^b_{iby}', '\it\omega^b_{ibz}');
msplot(122, tt(2:end), vm/ts, '\itf^b\rm_{sf} / ( m.s^{-2} )');
      legend('\itf^b\rm_{sf\itx}', '\itf^b\rm_{sf\ity}', '\itf^b\rm_{sf\itz}');
% ��ά�켣
%%
figure
delpos=deltapos(pos);
plot3(pos(:,1)/arcdeg,pos(:,2)/arcdeg,delpos(:,3));
% axis([0,-8000,]);
grid on;
xlabel('γ��/��'),ylabel('����/��'),zlabel('��/m');
%% ������Ϣ
wim0=my_getII(wm(1,:)',ts);
wim=[wim0,my_getII(wm',ts)];%wimע������������ϵ�µ�����
fm=[my_getII(vm(1,:)',ts),my_getII(vm',ts)];
MINS.wim=wim;
MINS.wim0=wim0;
MINS.fm=fm;
MINS.ts=ts;
fprintf('���߷�����ɣ�\n');
%% �ӹ��˶��������
fprintf('�ӹ߷���...'); 
Sinf.ts=MINS.ts;
len=size(wm,1)+1;
Sinf.len=size(wm,1)+1;% len=7200;
%�ӹ��˶�����
Sinf.Rlist=[2 4;
            0 0;
             0 0];
Sinf.ulist=[1 1]*arcdeg;
Sinf.flist=[1 1]*0.01;
Sinf.aerrlist=[1 1;
          2 2;
          3 3]*arcdeg;
Sinf.wim0=wim0;
SmoveCell=my_nSmovePack(Sinf);
Smove1=SmoveCell{1};
Smove2=SmoveCell{2};
% Smove1.um=1*arcdeg;%Rλ�ô���Ч�Ƕ����
% Smove1.f=0.01;%��Ƶ��
% Smove1.u=[zeros(1,len);Smove1.um*sin(2*pi*Smove1.f*(0:ts:((len-1)*ts)));zeros(1,len)];
% Smove1.R0=[2;0;0];%����ʱ ���Ӽ�0ʱ��ǰ��ʼ���λ��
% Smove1.Rf0=[0;0;0];%����ʱ ��0ʱ��ǰƫ��
% Smove1.Rf=[-Smove1.R0(1)*sin(Smove1.u(2,:)).*sin(Smove1.u(2,:));zeros(1,len);Smove1.R0(1)*sin(Smove1.u(2,:)).*cos(Smove1.u(2,:))];
% Smove1.R=Smove1.R0+Smove1.Rf;
% Smove1.att=2*Smove1.u;
% Smove1.att_s0=[0;0;0]*arcdeg;%����ʱ ���Ӽ�0ʱ��ǰ��ʼ��ԽǶ�
% Smove1.diffRf0=[0;0;0];%Rf�ĳ�ʼ�仯��
% Smove1.U0=Smove1.diffRf0+cross(wim0,Smove1.R0);
% Smove1.atterr0=[1;2;3]*arcdeg;%ע��������ת˳��zxy
% %�ӹ� 2 ���˶�
% Smove2=Smove1;
% Smove2.um=1*arcdeg;%Rλ�ô���Ч�Ƕ����
% Smove2.f=0.01;%��Ƶ��
% Smove2.u=[zeros(1,len);Smove2.um*sin(2*pi*Smove2.f*(0:ts:((len-1)*ts)));zeros(1,len)];
% Smove2.R0=[4;0;0];%����ʱ ���Ӽ�0ʱ��ǰ��ʼ���λ��
% Smove2.Rf0=[0;0;0];%����ʱ ��0ʱ��ǰƫ��
% Smove2.Rf=[-Smove2.R0(1)*sin(Smove2.u(2,:)).*sin(Smove2.u(2,:));zeros(1,len);Smove2.R0(1)*sin(Smove2.u(2,:)).*cos(Smove2.u(2,:))];
% Smove2.R=Smove2.R0+Smove2.Rf;
% Smove2.att=2*Smove2.u;
% Smove2.att_s0=[0;0;0]*arcdeg;%����ʱ ���Ӽ�0ʱ��ǰ��ʼ��ԽǶ�
% Smove2.diffRf0=[0;0;0];%Rf�ĳ�ʼ�仯��
% Smove2.U0=Smove2.diffRf0+cross(wim0,Smove2.R0);
% Smove2.atterr0=[1;2;3]*arcdeg;%ע��������ת˳��zxy
%% �ӹ߷���������
SINS1=my_invRIpack(MINS,Smove1);
SINS2=my_invRIpack(MINS,Smove2);
%% �Ƕ���ֵ zxyת��
fprintf('�ӹ߷�����ɣ�\n���ӹ���Ե���...%5.0f %%',0); 
%% �ӹ�������
SERR1.eb = [0.1;0.1;0.1]*dph; SERR1.web = [0.1;0.1;0.1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR1.db = [20;20;20]*ug; SERR1.wdb = [20;20;20]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
SERR2.eb = 0*[0.1;0.1;0.1]*dph; SERR2.web = 0*[0.1;0.1;0.1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR2.db = 0*[20;20;20]*ug; SERR2.wdb = 0*[20;20;20]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
%% ��Ե���
atterr0_1=[1.10;2.01;3.50]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
atterr0_2=[1.10;2.01;3.50]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
SINSR1=my_SINSgetResult(MINS,SINS1,Sinf.len/2,SERR1,atterr0_1);
SINSR2=my_SINSgetResult(MINS,SINS2,Sinf.len/2,SERR2,atterr0_2);
fprintf('������Ե�����ɣ�\n�˲�...%5.0f %%',0);
%% �˲�
nts = 2*ts; %�������Ͳ���ʱ��
%�ӹ�1
KFinit1.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit1.rk = [0.0001;0.0001;0.0001];  
KFinit1.Rk = diag(KFinit1.rk)^2;
KFinit1.P0 = diag([[0.1;0.1;10]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR1,~]=my_getFResult(MINS,SINS1,KFinit1,atterr0_1,Sinf.len/2,...
    SINSR1.ws_m_addnoise,SINSR1.fs_m_addnoise);
%�ӹ�2
KFinit2.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit2.rk = [0.0001;0.0001;0.0001];  
KFinit2.Rk = diag(KFinit2.rk)^2;
KFinit2.P0 = diag([[0.1;0.1;10]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR2,~]=my_getFResult(MINS,SINS2,KFinit2,atterr0_2,Sinf.len/2,...
    SINSR2.ws_m_addnoise,SINSR2.fs_m_addnoise);
fprintf('�˲���ɣ�\n'); 
%% 
R=SINS1.R;%����λ��
figure%�������λ���봿�������λ�öԱ�
subplot(311)
T=ts:ts:(Sinf.len*ts);
plot(T(2:4:end),R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.Rall(1:2:end,1),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('x/m');
title('���ӹߵ����λ��R');
grid on;
subplot(312)
plot(T(2:4:end),R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.Rall(1:2:end,2),'--','LineWidth',2)
grid on;
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('y/m');
title('���ӹߵ����λ��R');
subplot(313)
plot(T(2:4:end),R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.Rall(1:2:end,3),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('z/m');
title('���ӹߵ����λ��R');
grid on;
%%
figure%���۽Ƕ��봿��Խ���ǶȶԱ�
subplot(311)
plot(T(2:4:end),(SINS1.atttrue(1,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(1,1:2:end)/arcdeg,'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('������/��');
title('���ӹߵ������̬');

grid on;
subplot(312)
plot(T(2:4:end),(SINS1.atttrue(2,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(2,1:2:end)/arcdeg,'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('��ת��/��');
title('���ӹߵ������̬');

grid on;
subplot(313)
plot(T(2:4:end),(SINS1.atttrue(3,2:4:end))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end),SINSR1.attall(3,1:2:end)/arcdeg,'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('�����/��');
title('���ӹߵ������̬');

grid on;
%% 
figure%�������λ���봿�������λ�öԱ�
subplot(311)
T=ts:ts:(Sinf.len*ts);
plot(T(2:4:end),SINS1.R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,1),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('x/m');
title('���ӹߵ����λ��R���˲���');
grid on;
subplot(312)
plot(T(2:4:end),SINS1.R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,2),'--','LineWidth',2)
grid on;
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('y/m');
title('���ӹߵ����λ��R���˲���');
subplot(313)
plot(T(2:4:end),SINS1.R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),SINSFR1.Rall(1:2:end,3),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('z/m');
title('���ӹߵ����λ��R���˲���');
grid on;
%%
figure%���۽Ƕ��봿��Խ���ǶȶԱ�
subplot(311)
plot(T(2:4:end-1),(SINS1.atttrue(1,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(1,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('������/��');
title('���ӹߵ������̬���˲���');

grid on;
subplot(312)
plot(T(2:4:end-1),(SINS1.atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(2,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('��ת��/��');
title('���ӹߵ������̬���˲���');

grid on;
subplot(313)
plot(T(2:4:end-1),(SINS1.atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
hold on;
plot(T(2:4:end-1),SINSFR1.attall(3,1:2:end-1)/arcdeg,'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('�����/��');
title('���ӹߵ������̬���˲���');

grid on;
%% 
figure
subplot(311)
T=ts:ts:(Sinf.len*ts);
plot(T(2:4:end),SINSFR1.Rall(1:2:end,1)'-SINS1.R(1,2:4:end),'LineWidth',2)
legend('X');
xlabel('ʱ��/s');ylabel('x/m');
title('���ӹߵ����λ�����');
grid on;
subplot(312)
plot(T(2:4:end),SINSFR1.Rall(1:2:end,2)'-SINS1.R(2,2:4:end),'LineWidth',2)
grid on;
legend('Y');
xlabel('ʱ��/s');ylabel('y/m');
title('���ӹߵ����λ�����');
subplot(313)
plot(T(2:4:end),SINSFR1.Rall(1:2:end,3)'-SINS1.R(3,2:4:end),'LineWidth',2)
legend('Z');
xlabel('ʱ��/s');ylabel('z/m');
title('���ӹߵ����λ�����');
grid on;
%%
figure
subplot(311)
plot(T(2:4:end-1),(SINSFR1.attall(1,1:2:end-1)-SINS1.atttrue(1,2:4:end-1))/arcdeg,'LineWidth',2)
xlabel('ʱ��/s');ylabel('������/��');title('���ӹߵ������̬���');grid on;
subplot(312)
plot(T(2:4:end-1),(SINSFR1.attall(2,1:2:end-1)-SINS1.atttrue(2,2:4:end-1))/arcdeg,'LineWidth',2)
xlabel('ʱ��/s');ylabel('��ת��/��');title('���ӹߵ������̬���');grid on;
subplot(313)
plot(T(2:4:end-1),(SINSFR1.attall(3,1:2:end-1)-SINS1.atttrue(3,2:4:end-1))/arcdeg,'LineWidth',2)
xlabel('ʱ��/s');ylabel('�����/��');title('���ӹߵ������̬���');grid on;
%% ͳ�����
errlen=1/(Smove1.f*nts);
Rerr=my_getSErr(SINSFR1.Rall(1:end,:)'-SINS1.R(:,1:2:end),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSFR1.attall(:,1:end)-SINS1.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('�˲�ͳ��,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSFR1.Rall(1:end-1,:)'-SINS1.R(:,1:2:end-2),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSFR1.attall(:,1:end-1)-SINS1.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('�˲���ֵ,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1/(Smove1.f*nts);
Rerr=my_getSErr(SINSR1.Rall(1:end,:)'-SINS1.R(:,1:2:end),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSR1.attall(:,1:end)-SINS1.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('����ͳ��,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSR1.Rall(1:end-1,:)'-SINS1.R(:,1:2:end-2),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSR1.attall(:,1:end-1)-SINS1.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('������ֵ,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
%%
disp "SINS2"
errlen=1/(Smove2.f*nts);
Rerr=my_getSErr(SINSFR2.Rall(1:end,:)'-SINS2.R(:,1:2:end),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSFR2.attall(:,1:end)-SINS2.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('�˲�ͳ��,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSFR2.Rall(1:end-1,:)'-SINS2.R(:,1:2:end-2),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSFR2.attall(:,1:end-1)-SINS2.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('�˲���ֵ,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1/(Smove2.f*nts);
Rerr=my_getSErr(SINSR2.Rall(1:end,:)'-SINS2.R(:,1:2:end),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSR2.attall(:,1:end)-SINS2.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('����ͳ��,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSR2.Rall(1:end-1,:)'-SINS2.R(:,1:2:end-2),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSR2.attall(:,1:end-1)-SINS2.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('������ֵ,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));