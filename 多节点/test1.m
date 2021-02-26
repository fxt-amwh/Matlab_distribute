%% һ��+���ӣ��߶�ģ�黯��ʹ�ýǶȷֱ�������Ķ�ڵ��ں�
% 1��������ƹ켣�˶����������ߴ��������
% 2��������ƻ����˶��������ӹߴ���������
% 3���ӹ�����ע�����
% 4����������Ե���
% 5���˲�
% 6��������
% �汾ʱ�䣺2020/12/19 ��ڵ��������
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
%         0,         0,          0,         0,         120                   %360
%         60,         0,          0,         0,         30                   %390
%         0,         0,          0,         0,         220                   %610
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
wat(:,1:3) = wat(:,1:3)*arcdeg/60;  % deg/min->deg/s
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
MINS.wim=wim;MINS.wim0=wim0;MINS.fm=fm;MINS.ts=ts;
fprintf('���߷�����ɣ�\n');
%% ȫ�����
flag.EnReli=false;%�Ƿ񴿽�����Ե������˲�
flag.EnFusion=true;%�Ƿ�����ں�
nameList=["SINS1","SINS2","SINS3","SINS4"];
% �ӹ��˶��������
fprintf('�ӹ߷���...'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;
Sinf.Rlist=[1 2 3 4;
            0 0 0 0;
            0 0 0 0];
Sinf.ulist=[1 2 3 4]*arcdeg;%�ӹߴ�������
Sinf.flist=[1 1 1 1]*0.01;%����Ƶ��
% Sinf.aerrlist=[ 1 1 1 1;
%                 2 2 2 2;
%                 3 3 3 3]*arcdeg;%�ӹ߰�װ����
Sinf.aerrlist=[ 0 0 0 0;
                0 0 0 0;
                0 0 0 0]*arcdeg;%�ӹ߰�װ����
Sinf.wim0=MINS.wim0;SinfCell{1,1}=Sinf;
Smove=my_nSmovePackN(SinfCell);
[M,N]=size(Smove);%��������ƶ���������
% �ӹ�������
SERRList.eb=zeros(3,N,M);%zeros(3,N,M)
SERRList.web=zeros(3,N,M);%zeros(3,N,M)
SERRList.db=zeros(3,N,M);%zeros(3,N,M)
SERRList.wdb=zeros(3,N,M);%zeros(3,N,M)
% for m=1:M
%     SERRList.eb(:,:,m)=[0.1 0.1 0.1;%
%                         0.1 0.1 0.1;%
%                         0.1 0.1 0.1]*dph;%���ݳ�ֵ��ƫ һ�ж�Ӧһ��
%     SERRList.web(:,:,m)=[0.1 0.1 0.1;%
%                          0.1 0.1 0.1;%
%                          0.1 0.1 0.1]*dpsh;%�Ƕ�������� һ�ж�Ӧһ��
%     SERRList.db(:,:,m)=[20 20 20;%
%                         20 20 20;%
%                         20 20 20]*ug;%���ٶȼƳ�ֵƫֵ һ�ж�Ӧһ��
%     SERRList.wdb(:,:,m)=[20 20 20;%
%                          20 20 20;%
%                          20 20 20]*ugpsHz;%�ٶ�������� һ�ж�Ӧһ��
% end
for m=1:M
    for n=1:N
        SERRList.eb(:,n,m)=[1;1;1]*dph;%���ݳ�ֵ��ƫ һ�ж�Ӧһ��
        SERRList.web(:,n,m)=[1;1;1]*dpsh;%�Ƕ�������� һ�ж�Ӧһ��
        SERRList.db(:,n,m)=[200;200;200]*ug;%���ٶȼƳ�ֵƫֵ һ�ж�Ӧһ��
        SERRList.wdb(:,n,m)=[200;200;200]*ugpsHz;%�ٶ�������� һ�ж�Ӧһ��
    end
end
% ��ʼʱ�İ�װ����Լֵ
atterrList=zeros(3,N,M);%zeros(3,N,M);
for m=1:M
%     atterrList(:,:,m)=[ 1.01 1.01 1.01 1.01;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
%                        2.01 2.01 2.01 2.01;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
%                        3.01 3.01 3.01 3.01]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
    atterrList(:,:,m)=[ 0.01 0.01 0.01 0.01;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
                        0.01 0.01 0.01 0.01;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
                        0.01 0.01 0.01 0.01]*arcdeg;%ע��������ת˳��zxy ��׼ȷ�ĳ�ʼ��
%     for n=1:N
%         Smove{m,n}.atterrList=[0;0;0];
%     end
                   
end
%% �ӹ߷���������
SINS=my_invRIpackN(MINS,Smove);
fprintf('�ӹ߷�����ɣ�\n'); 
%% �ӹ����ת��ΪԪ��
SERR=my_getSERRN(SERRList);
%% ��Ե���
atterr=my_getAtterr(atterrList);
if flag.EnReli
    SINSR=my_SINSgetResultN(MINS,SINS,SinfCell{1,1}.len/2,SERR,atterr); 
end
%% �˲�
nts = 2*ts; %�������Ͳ���ʱ��
KFinit=my_KFinitN(SERR,SINS,2*ts);
if flag.EnReli
    % [SINSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,1);%�����ӹ߷����
    flag.select=1;
    if flag.EnFusion
        [SINSFR,Filter,NSINSFR]=my_getFResultLoopN_new(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,flag);%���м��а�
    else
        [SINSFR,Filter,~]=my_getFResultLoopN_new(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SINSR,flag);%���м��а�
    end
else
    flag.select=2;
    % [SINSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,2);%�����ӹ߷����
    if flag.EnFusion
        [SINSFR,Filter,NSINSFR]=my_getFResultLoopN_new(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,flag);%���м��а�
    else
        [SINSFR,Filter,~]=my_getFResultLoopN_new(MINS,SINS,KFinit,atterr,SinfCell{1,1}.len/2,SERR,flag);%���м��а�
    end
end
%% ������Ե����ӹ�1λ�ý��
flag.Ceil=1;flag.flagSINS_M=1;flag.flagSINS_N=1;flag.figureFlag=1;myfigure;
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
if flag.EnFusion
    % �����ӹ�1�˲�λ�����
    flag.figureFlag=7;myfigure;
    % �����ӹ�1�˲���̬���
    flag.figureFlag=8;myfigure;
end

%% ͳ�����
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
if flag.EnFusion
    disp('�ںϽ��');
    errlen=1/(Smove{1,1}.f*nts);
    Rerr=my_getSErr(NSINSFR.Rall(1:end,:)'-SINS{1,1}.R(:,1:2:end),errlen,2)*1000;        %λ����λmm
    Atterr=my_getSErr(NSINSFR.attnewall(:,1:end)-SINS{1,1}.attuftrue(:,1:2:end),errlen,2)/arcdeg*60; %��̬��λ��
    fprintf('�˲�ͳ��,λ�����%f %f %f��mm������̬���%f %f %f���֣�\n',...
        Rerr(1),Rerr(2),Rerr(3),Atterr(1),Atterr(2),Atterr(3));
    errlen=1;
    beg=ceil(200/nts);
    End=ceil(500/nts);
    Rerr=my_getSErr(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),errlen,2)*1000;        %λ����λmm
    Rserr=my_getVErr(NSINSFR.Rall(1:end-1,:)'-SINS{1,1}.R(:,1:2:end-2),beg,End,2)*1000;        %λ����λmm
    Atterr=my_getSErr(NSINSFR.attnewall(:,1:end-1)-SINS{1,1}.attuftrue(:,1:2:end-2),errlen,2)/arcdeg*60; %��̬��λ��
    Attserr=my_getVErr(NSINSFR.attnewall(:,1:end-1)-SINS{1,1}.attuftrue(:,1:2:end-2),beg,End,2)/arcdeg*60;
    fprintf('�˲���ֵ,λ�����%f %f %f��mm��,λ������%f %f %f(mm)����̬���%f %f %f���֣�,��̬������%f %f %f(��)\n',...
        Rerr(1),Rerr(2),Rerr(3),Rserr(1),Rserr(2),Rserr(3),Atterr(1),Atterr(2),Atterr(3),Attserr(1),Attserr(2),Attserr(3));   
end