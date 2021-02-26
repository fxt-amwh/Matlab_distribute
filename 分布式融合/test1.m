%% һ��+���� ���Խ������ӹ�����ת������һ�ӹߴ�
% ���ݵ�ת��ֻ��Ҫ��ת
% �ӼƵ�ת������Ҫ���Ǹ˱ۼ��ٶ�
% �汾ʱ�䣺2021/02/06 һ�����Ӷ�ڵ��ںϳ�������
clc
clear
close all
gvar;    % ����ȫ�ֱ���
ts = 0.05;%����ʱ��
g2d=180/pi;
d2g=1/g2d;
att0 = [0;0;0]*arcdeg; vn0 = [0;200;0]; pos0 = [[34;108]*arcdeg;105];
%     ���������� ��������� ��λ������ ������ٶ� ����ʱ�� �����ι켣��ƣ�
wat = [  0,         0,          0,         0,         10       %��ֹ
        60,         0,          0,         0,         60       %����
        -60,        0,          0,         0,         60       %����
        0,          0,          0,         0,         50       %����
        0,          0,          60,       0,         360      %����ת��
        0,          0,          0,         0,         50       %����
        -60,        0,          0,         0,         60       %����
        60,         0,          0,         0,        60       %����
        0,          0,          0,         0,         10       %����
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
%% ��Զ�ӹ���Ϣ
simf=0.01;
uf_lw_m=5*d2g;
R_lw=[2;0;0];
%% ȫ�����
flag.EnReli=false;%�Ƿ񴿽�����Ե������˲�
flag.EnFusion=true;%�Ƿ�����ں�
nameList=["SINS1_1","SINS1_2","SINS1_3","SINS1_4"
    "SINS2_1","SINS2_2","SINS2_3","SINS2_4"];
% �ӹ��˶��������
fprintf('�ӹ߷���...\n'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;


% Sinf.Rlist=[2 2 2 2 2;
%             0 0 0 0 0;
%             0 0 0 0 0];
Sinf.Rlist=[linspace(1,R_lw(1,1),4);zeros(2,4)];
Sinf.ulist=uf_lw_m*Sinf.Rlist(1,:)/R_lw(1,1);%�ӹߴ�������
Sinf.flist=simf*ones(1,size(Sinf.Rlist,2));%����Ƶ��
Sinf.aerrlist=[ 1  1  1  1 ;
                2  2  2  2 ;
                3  3  3  3 ]*d2g;%�ӹ߰�װ����
Sinf.wim0=MINS.wim0;SinfCell{1,1}=Sinf;
% Sinf.Rlist(2,:)=0.5;%�ڶ���y��λ��
SinfCell{1,2}=Sinf;
Smove=my_nSmovePackN(SinfCell);
[M,N]=size(Smove);%��������ƶ���������
% �ӹ�������
SERRList.eb=zeros(3,N,M);%zeros(3,N,M)
SERRList.web=zeros(3,N,M);%zeros(3,N,M)
SERRList.db=zeros(3,N,M);%zeros(3,N,M)
SERRList.wdb=zeros(3,N,M);%zeros(3,N,M)
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
    atterrList(:,:,m)=[ 2  1  2  1 ;
                        2  3  2  2 ;
                        3  3  3  2 ]*d2g;%�ӹ߰�װ����
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