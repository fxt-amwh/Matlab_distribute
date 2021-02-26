clc;
clear;
close all;
%% ������Ϣ
g2d=180/pi;
d2g=1/g2d;
ts=0.01;
simT=10;
simf=2;
uf_lw_m=5*d2g;
R_lw=[5;0;0];
t=0:ts:(simT/simf);
wim0=zeros(3,1);
wim=[wim0,zeros(3,length(t)-1)];%wimע������������ϵ�µ�����
fm=wim;
MINS.wim=wim;MINS.wim0=wim0;MINS.fm=fm;MINS.ts=ts;
fprintf('���߷�����ɣ�\n');
%% ȫ�����
nameList=["SINS1","SINS2","SINS3","SINS4"];
% �ӹ��˶��������
fprintf('�ӹ߷���...\n'); 
Sinf.ts=MINS.ts;
len=size(MINS.wim,2);
Sinf.len=size(MINS.wim,2);% len=7200;


% Sinf.Rlist=[2 2 2 2;
%             0 0 0 0;
%             0 0 0 0];
Sinf.Rlist=[linspace(1,R_lw(1,1),5);zeros(2,5)];
Sinf.ulist=uf_lw_m*Sinf.Rlist(1,:)/R_lw(1,1);%�ӹߴ�������
Sinf.flist=simf*ones(1,size(Sinf.Rlist,2));%����Ƶ��
Sinf.aerrlist=[ 1 1 1 1 1;
                2 2 2 2 2;
                3 3 3 3 3]*d2g;%�ӹ߰�װ����
Sinf.wim0=MINS.wim0;SinfCell{1,1}=Sinf;
Sinf.Rlist(2,:)=1;
SinfCell{1,2}=Sinf;
Smove=my_nSmovePackN(SinfCell);
[M,N]=size(Smove);%��������ƶ���������
figure
now_t=1;
x1=zeros(M,N+1);
p1=zeros(M,N+1);
y1=zeros(M,N+1);
for i=1:len
    for m=1:M
        for n=1:N
            x1(m,n+1)=Smove{m,n}.R(1,now_t);
            y1(m,n+1)=Smove{m,n}.R(2,now_t);
            p1(m,n+1)=Smove{m,n}.R(3,now_t); 
        end
    end
%     plot(x1(1,:),p1(1,:),'-*');
%     y1(:,2:end)=[0*ones(1,N);ones(1,N)];��
    maxx=max(R_lw(1,:));
    if M>1
        surf(x1,y1,p1);
        title("�����ӹ������˶�����")
        hold on;
        plot3(x1,y1,p1,"*");
        hold off;
        text(-0.5,0,0,"����")
        axis([0,maxx,-maxx/2,maxx/2,-maxx/2,maxx/2])
    else
        plot(x1(1,:),p1(1,:),'-*');
        title("�����ӹ������˶�����");
        text(-0.5,0,"����")
        axis([0,maxx,-maxx/2,maxx/2])
    end
    pause(ts);
    now_t=now_t+1;
end
