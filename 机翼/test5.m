%% һ��+���ӣ��߶�ģ�黯-�˲��赥�����ã��㱨����2
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
ts = 0.01;%����ʱ��
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;105];
%     ���������� ��������� ��λ������ ������ٶ� ����ʱ�� �����ι켣��ƣ�
% wat = [  0,         -60,          120,         0,         30               %30
%         0,         60,          -120,        0,         60                 %90
%         0,         0,          0,         0,         120                   %210
%         -60,         0,          0,         0,         30                  %240
% %         0,         0,          360*60,         0,         360       %����ת��
%         0,         0,          0,         0,         120                   %360
%         60,         0,          0,         0,         30                   %390
%         0,         0,          0,         0,         220                   %610
%         ];    %��ֹ
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
flag.EnReli=false;%�Ƿ񴿽�����Ե������˲�
flag.EnFusion=true;%�Ƿ�����ں�
nameList=["SINS1","SINS2","SINS3"];
fprintf('�ӹ߷���...'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;
Sinf.Rlist=[5 5 5;
            0 0 0;
            0 0 0];
Sinf.ulist=[1 2 3]*arcdeg;
Sinf.flist=[1 1 1]*0.01;
Sinf.aerrlist=[ 1 1 1;
                2 2 2;
                3 3 3]*arcdeg;
Sinf.wim0=MINS.wim0;
SinfCell{1,1}=Sinf;
% SmoveCell=my_nSmovePack(Sinf);
Smove=my_nSmovePackN(SinfCell);
%% �ӹ߷���������
SINS=my_invRIpackN(MINS,Smove);
%% �Ƕ���ֵ zxyת��
fprintf('�ӹ߷�����ɣ�\n'); 
%% �ӹ�������
SERR1.eb = [1;1;1]*dph; SERR1.web = [1;1;1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR1.db = [200;200;200]*ug; SERR1.wdb = [200;200;200]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
SERR2.eb = [1;1;1]*dph; SERR2.web = [1;1;1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR2.db = [200;200;200]*ug; SERR2.wdb = [200;200;200]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
SERR3.eb = [1;1;1]*dph; SERR3.web = [1;1;1]*dpsh;   %���ݳ�ֵ��ƫ���Ƕ��������
SERR3.db = [200;200;200]*ug; SERR3.wdb = [200;200;200]*ugpsHz;  %���ٶȼƳ�ֵƫֵ���ٶ��������
SERR{1,1}=SERR1;SERR{1,2}=SERR2;SERR{1,3}=SERR3;
%% ��Ե���
atterr0_1=[1.01;2.01;3.01]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
atterr0_2=[1.01;2.01;3.01]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
atterr0_3=[1.01;2.01;3.01]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
atterr{1,1}=atterr0_1;atterr{1,2}=atterr0_2;atterr{1,3}=atterr0_3;
SINSR=my_SINSgetResultN(MINS,SINS,SinfCell{1,1}.len/2,SERR,atterr);
%% �˲�
nts = 2*ts; %�������Ͳ���ʱ��
%�ӹ�1
KFinit1.Qk = diag([SERR{1,1}.web; SERR{1,1}.wdb;])^2*nts;
KFinit1.rk = [0.001;0.001;0.001]*sqrt(SINS{1,1}.R0'*SINS{1,1}.R0);  
KFinit1.Rk = diag(KFinit1.rk)^2;
KFinit1.P0 = diag([[0.1;0.1;0.1]*arcdeg; [20;20;20]; [0.002;0.002;0.002]*sqrt(SINS{1,1}.R0'*SINS{1,1}.R0);
         SERR{1,1}.eb; SERR{1,1}.db])^2;
KFinit{1,1}=KFinit1;

%�ӹ�2
KFinit2.Qk = diag([SERR{1,2}.web; SERR{1,2}.wdb;])^2*nts;
KFinit2.rk = [0.001;0.001;0.001]*sqrt(SINS{1,2}.R0'*SINS{1,2}.R0);   
KFinit2.Rk = diag(KFinit2.rk)^2;
KFinit2.P0 = diag([[0.1;0.1;0.1]*arcdeg; [20;20;20]; [0.002;0.002;0.002]*sqrt(SINS{1,2}.R0'*SINS{1,2}.R0);
         SERR{1,2}.eb; SERR{1,2}.db])^2;
KFinit{1,2}=KFinit1;
%�ӹ�3
KFinit3.Qk = diag([SERR{1,3}.web; SERR{1,3}.wdb;])^2*nts;
KFinit3.rk = [0.001;0.001;0.001]*sqrt(SINS{1,3}.R0'*SINS{1,3}.R0);  
KFinit3.Rk = diag(KFinit3.rk)^2;
KFinit3.P0 = diag([[0.1;0.1;0.1]*arcdeg; [20;20;20]; [0.002;0.002;0.002]*sqrt(SINS{1,3}.R0'*SINS{1,3}.R0);
         SERR{1,3}.eb; SERR{1,3}.db])^2;
KFinit{1,3}=KFinit1;
[SINSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,1);
%% ������Ե����ӹ�1λ�ý��
flag.EnReli=true;
flag.Ceil=1;flag.flagSINS_M=1;flag.flagSINS_N=3;flag.figureFlag=1;myfigure;
% ������Ե����ӹ�1��̬���
flag.figureFlag=2;myfigure;
% �����ӹ�1�˲�λ�ý��
flag.figureFlag=3;myfigure;
% �����ӹ�1�˲���̬���
flag.figureFlag=4;myfigure;
% �����ӹ�1�˲�λ�����
flag.figureFlag=5;myfigure;
% �����ӹ�1�˲���̬���
flag.figureFlag=6;myfigure;
%% ͳ�����
M=1;N=3;
for m=1:M
    for n=1:N
        disp(nameList(m,n));
        errlen=1/(Smove{m,n}.f*nts);
        Rerr=my_getSErr(SINSFR{m,n}.Rall(1:end,:)'-SINS{m,n}.R(:,1:2:end),errlen,2)*1000;        %λ����λmm
        Atterr=my_getSErr(SINSFR{m,n}.attall(:,1:end)-SINS{m,n}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %��̬��λ��
        fprintf('�˲�ͳ��,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
            Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
        errlen=1;
        beg=ceil(200/nts);
        End=ceil(500/nts);
        Rerr=my_getSErr(SINSFR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),errlen,2)*1000;        %λ����λmm
        Rserr=my_getVErr(SINSFR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),beg,End,2)*1000;        %λ����λmm
        Atterr=my_getSErr(SINSFR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %��̬��λ��
        Attserr=my_getVErr(SINSFR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),beg,End,2)/arcdeg*60;
        fprintf('�˲���ֵ,λ�����%f %f %f��mm��,λ������%f %f %f(mm)����̬���%f %f %f���֣�,��̬������%f %f %f(��)\n',...
            Rerr(1),Rerr(2),Rerr(3),Rserr(1),Rserr(2),Rserr(3),Atterr(1),Atterr(2),Atterr(3),Attserr(1),Attserr(2),Attserr(3));
        if flag.EnReli
            errlen=1/(Smove{m,n}.f*nts);
            Rerr=my_getSErr(SINSR{m,n}.Rall(1:end,:)'-SINS{m,n}.R(:,1:2:end),errlen,2)*1000;        %λ����λmm
            Atterr=my_getSErr(SINSR{m,n}.attall(:,1:end)-SINS{m,n}.atttrue(:,1:2:end),errlen,2)/arcdeg*60; %��̬��λ��
            fprintf('����ͳ��,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
                Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
            errlen=1;
            Rerr=my_getSErr(SINSR{m,n}.Rall(1:end-1,:)'-SINS{m,n}.R(:,1:2:end-2),errlen,2)*1000;        %λ����λmm
            Atterr=my_getSErr(SINSR{m,n}.attall(:,1:end-1)-SINS{m,n}.atttrue(:,1:2:end-2),errlen,2)/arcdeg*60; %��̬��λ��
            fprintf('������ֵ,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
                Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
        end
    end
end
