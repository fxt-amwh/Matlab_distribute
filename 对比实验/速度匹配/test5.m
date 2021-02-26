%�ļ����ƣ�test5
%������Ϣ��������
%��������������������ݽṹ���Լ��ӹ����ݽṹ�� ������ ����Ϊ�˷�����Ը����ڲ�����
%�汾ʱ�䣺2021/1/30 19:17
%% �����˶�������Ʒ���
clc
clear
close all
gvar;    % ����ȫ�ֱ���
ts = 0.01;%����ʱ��
g2d=180/pi;
d2g=1/g2d;
att0 = [0;0;0]*arcdeg; vn0 = [0;100;0]; pos0 = [[34;108]*arcdeg;100];
%     ���������� ��������� ��λ������ ������ٶ� ����ʱ�� �����ι켣��ƣ�
% wat = [  0,         -60,          120,         0,         30               %30
%         0,         60,          -120,        0,         60                 %90
%         0,         0,          0,         0,         120                   %210
%         -60,         0,          0,         0,         30                  %240
%         0,         0,          0,         0,         120                   %360
%         60,         0,          0,         0,         30                   %390
%         0,         0,          0,         0,         220                   %610
%         ];    %��ֹ
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
% ���TR
TR.att0=att0;
TR.vn0=vn0;
TR.pos0=pos0;
TR.vn=vn;
TR.attm=attm;
TR.pos=pos;
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
nameList=["SINS1_1","SINS1_2"];
[M,N]=size(nameList);%��������ƶ���������
% �ӹ��˶��������
fprintf('�ӹ߷���...\n'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;


Sinf.Rlist=[R_lw(1,1)*ones(1,N);zeros(2,N)];
Sinf.ulist=uf_lw_m*Sinf.Rlist(1,:)/R_lw(1,1);%�ӹߴ�������
Sinf.ulist(:,:)=0;%�ӹߴ�������
Sinf.flist=simf*ones(1,size(Sinf.Rlist,2));%����Ƶ��
Sinf.aerrlist=0*[ 1  1 ;
                2  2 ;
                3  3 ]*d2g;%�ӹ߰�װ����
Sinf.wim0=MINS.wim0;SinfCell{1,1}=Sinf;
 
Smove=my_nSmovePackN(SinfCell);
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
%         SERRList.eb(:,n,m)=[0.1;0.1;0.1]*dph;%���ݳ�ֵ��ƫ һ�ж�Ӧһ��
%         SERRList.web(:,n,m)=[0.1;0.1;0.1]*dpsh;%�Ƕ�������� һ�ж�Ӧһ��
%         SERRList.db(:,n,m)=[20;20;20]*ug;%���ٶȼƳ�ֵƫֵ һ�ж�Ӧһ��
%         SERRList.wdb(:,n,m)=[20;20;20]*ugpsHz;%�ٶ�������� һ�ж�Ӧһ��
%         SERRList.eb(:,n,m)=[0.5;0.5;0.5]*dph;%���ݳ�ֵ��ƫ һ�ж�Ӧһ��
%         SERRList.web(:,n,m)=[0.5;0.5;0.5]*dpsh;%�Ƕ�������� һ�ж�Ӧһ��
%         SERRList.db(:,n,m)=[100;100;100]*ug;%���ٶȼƳ�ֵƫֵ һ�ж�Ӧһ��
%         SERRList.wdb(:,n,m)=[100;100;100]*ugpsHz;%�ٶ�������� һ�ж�Ӧһ��
        SERRList.eb(:,n,m)=[1;1;1]*dph;%���ݳ�ֵ��ƫ һ�ж�Ӧһ��
        SERRList.web(:,n,m)=[1;1;1]*dpsh;%�Ƕ�������� һ�ж�Ӧһ��
        SERRList.db(:,n,m)=[200;200;200]*ug;%���ٶȼƳ�ֵƫֵ һ�ж�Ӧһ��
        SERRList.wdb(:,n,m)=[200;200;200]*ugpsHz;%�ٶ�������� һ�ж�Ӧһ��
    end
end
% ��ʼʱ�İ�װ����Լֵ
atterrList=zeros(3,N,M);%zeros(3,N,M);
for m=1:M
    atterrList(:,:,m)=0*[ 1  1 ;
                        2  2 ;
                        3  3 ]*d2g;%�ӹ߰�װ����
end
%% �ӹ߷���������
SINS=my_invRIpackN(MINS,Smove);
fprintf('�ӹ߷�����ɣ�\n'); 
%% �ӹ����ת��ΪԪ��
SERR=my_getSERRN(SERRList);
%% �ӹ����ע��
len=length(SINS{1,2}.wis);
for i=1:len/2
    [ws_m, fs_m] = my_imuadderr(SINS{1,2}.wis(:,(2*i-1):(2*i))', SINS{1,2}.fs(:,(2*i-1):(2*i))', SERR{1,2}.eb, SERR{1,2}.web, SERR{1,2}.db, SERR{1,2}.wdb, ts);%�ӹ�ע������
end
%%
% qnb=a2qua(att0);
% vn=vn0;
% posSINS=pos0;
% % posSINS(3)=110;
% len=length(vm)-1;
% nn=2;nts=nn*ts;
% SINS=zeros(floor(len/2),10);kk=1;t=0;
% % eb=[0.01;0.01;0.01]*dph;
% % web=[0.001;0.001;0.001]*dpsh;
% % db=[10;10;10]*ug;
% % wdb=[0;0;0]*ugpsHz;
% eb=[0;0;0]*dph;
% web=[0;0;0]*dpsh;
% db=[0;0;0]*ug;
% wdb=[0;0;0]*ugpsHz;
% sensor1=zeros(2,6,len/2);
% sensor=zeros(2,6,len/2);
% num=1;
% for k=1:nn:len
%     t=t+nts;
%     [wm1,vm1]=imuadderr(wm(k:k+1,:),vm(k:k+1,:), eb, web, db, wdb, ts);
%     sensor1(:,:,num)=[wm1,vm1];sensor(:,:,num)=[wm1,vm1];
%     num=num+1;
%     [qnb,vn,posSINS]=my_insupdate(qnb,vn,posSINS,wm1,vm1,ts);
%     SINS(kk,:)=[q2att(qnb);vn;posSINS;t];kk=kk+1;
%     
%     if mod(t,100)<nts,disp(fix(t));end
% end
% %
% figure
% delpos=deltapos(pos);
% plot3(delpos(:,1),delpos(:,2),delpos(:,3),'LineWidth',5);
% hold on;
% delposSINS=deltapos(SINS(:,7:9));
% plot3(delposSINS(:,1),delposSINS(:,2),delposSINS(:,3));
% grid on;
% xlabel('γ�Ⱦ���/m'),ylabel('���Ⱦ���/m'),zlabel('��/m');
% err=delposSINS(end,:)-delpos(end,:);
% disp 'γ����� ������� �߶����'
% disp(err)