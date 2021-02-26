%% һ��+2�ӣ��м䱸�ݣ�
% 1��������ƹ켣�˶����������ߴ��������
% 2��������ƻ����˶��������ӹߴ���������
% 3���ӹ�����ע�����
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
wat = [  0,         -60,          120,         0,         30               %30
        0,         60,          -120,        0,         60                 %90
        0,         0,          0,         0,         120                   %210
        -60,         0,          0,         0,         30                  %240
%         0,         0,          360*60,         0,         360       %����ת��
        0,         0,          0,         0,         120                   %360
        60,         0,          0,         0,         30                   %390
        0,         0,          0,         0,         220                   %610
        
        0,         -60,          120,         0,         30                %640
        0,         60,          -120,        0,         60                 %700
        0,         0,          0,         0,         120                   %820
        -60,         0,          0,         0,         30                  %850
%         0,         0,          360*60,         0,         360            %����ת��
        0,         0,          0,         0,         120                   %970
        60,         0,          0,         0,         30                   %1000
        0,         0,          0,         0,         220                   %1220
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
figure
delpos=deltapos(pos);
plot3(delpos(:,1),delpos(:,2),delpos(:,3));
% axis([0,-8000,]);
grid on;
xlabel('γ�Ⱦ���/m'),ylabel('���Ⱦ���/m'),zlabel('��/m');
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
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;
Sinf.Rlist=[2 2 2;
            0 0 0;
            0 0 0];
Sinf.ulist=[1 2.5 2.5]*arcdeg;
Sinf.flist=[1 1 1]*0.01;
Sinf.aerrlist=0*[1 1 1;
          2 2 2;
          3 3 3]*arcdeg;
Sinf.wim0=MINS.wim0;
SmoveCell=my_nSmovePack(Sinf);
Smove1=SmoveCell{1};
Smove2=SmoveCell{2};
Smove3=SmoveCell{3};
%% �ӹ߷���������
SINS1=my_invRIpack(MINS,Smove1);
SINS2=my_invRIpack(MINS,Smove2);
SINS3=my_invRIpack(MINS,Smove3);
%% �Ƕ���ֵ zxyת��
fprintf('�ӹ߷�����ɣ�\n�˲�...%5.0f %%',0);
%% �ӹ�������
SERR1.eb = [0.1;0.1;0.1]*dph; SERR1.web = [0.1;0.1;0.1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR1.db = [20;20;20]*ug; SERR1.wdb = [20;20;20]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
SERR2.eb = 0*[0.1;0.1;0.1]*dph; SERR2.web = 0*[0.1;0.1;0.1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR2.db = 0*[20;20;20]*ug; SERR2.wdb = 0*[20;20;20]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
SERR3.eb = 0*[0.1;0.1;0.1]*dph; SERR3.web = 0*[0.1;0.1;0.1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR3.db = 0*[20;20;20]*ug; SERR3.wdb = 0*[20;20;20]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
%% �˲�
atterr0_1=0.1*[1.00;2.00;3.00]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
atterr0_2=0.1*[1.00;2.00;3.00]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
atterr0_3=0.1*[1.00;2.00;3.00]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
nts = 2*ts; %�������Ͳ���ʱ��
%�ӹ�1
KFinit1.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit1.rk = [0.0001;0.0001;0.0001];  
KFinit1.Rk = diag(KFinit1.rk)^2;
KFinit1.P0 = diag([[0.1;0.1;10]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR1,Filter1]=my_getFResult(MINS,SINS1,KFinit1,atterr0_1,Sinf.len/2,...
    SERR1);
%�ӹ�2
KFinit2.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit2.rk = [0.0001;0.0001;0.0001];  
KFinit2.Rk = diag(KFinit2.rk)^2;
KFinit2.P0 = diag([[0.1;0.1;10]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR2,Filter2]=my_getFResult(MINS,SINS2,KFinit2,atterr0_2,Sinf.len/2,...
    SERR2);
%�ӹ�3
KFinit3.Qk = diag([SERR1.web; SERR1.wdb;])^2*nts;
KFinit3.rk = [0.0001;0.0001;0.0001];  
KFinit3.Rk = diag(KFinit3.rk)^2;
KFinit3.P0 = diag([[0.1;0.1;10]*arcdeg; [20;20;20]; [0.01;0.01;0.01];
         [0.1;0.1;0.1]*dph; [20;20;20]*ug])^2;
[SINSFR3,Filter3]=my_getFResult(MINS,SINS3,KFinit3,atterr0_3,Sinf.len/2,...
    SERR3);
fprintf('�˲���ɣ�\n'); 
%% 
R=SINS1.R;%����λ��
%% 
figure%�������λ�����˲��������λ�öԱ�
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
figure%���۽Ƕ����˲�����ǶȶԱ�
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
%% λ�����
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
%% ��̬���
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
%%
disp "SINS3"
errlen=1/(Smove3.f*nts);
Rerr=my_getSErr(SINSFR3.Rall(1:end,:)'-SINS3.R(:,1:2:end),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSFR3.attall(:,1:end)-SINS3.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('�˲�ͳ��,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
errlen=1;
Rerr=my_getSErr(SINSFR3.Rall(1:end-1,:)'-SINS3.R(:,1:2:end-2),errlen,2)*1000;        %λ����λmm
Atterr=my_getSErr(SINSFR3.attall(:,1:end-1)-SINS3.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %��̬��λ��
fprintf('�˲���ֵ,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
    Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));