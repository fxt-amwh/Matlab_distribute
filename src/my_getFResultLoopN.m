function [INSSFR,Filter,NINSS]=my_getFResultLoopN(MINS,SINS,KFinit,atterr,len,selcCell,flag)
% �������ƣ�my_getFResultLoopN 
% �������ܣ����ӹ����������Ե����˲�
% ���룺MINS      :����
%       SINS      :�ӹ�
%       KFinit    :�ӹ��˲�����ʼֵ
%       atterr    :�ӹ߰�װ���ǹ��Ƴ�ʼֵ
%       len       :�˲�����       
%       selcCell  :select=0��selcCell���� �ӿڻ�����������
%                  select=1��selcCell��ʾ SINSR �ӹߴ���Ե������
%                  select=2��selcCell��ʾ SERR  �ӹ�����ʼֵ
%       select    :ѡ���ӹ������Դ
% �����INSSFR:ÿ���ӹ߽��Ԫ��
%       Filter:ÿ���ӹ��˲����м�ֵԪ��
%       NINSS :�ںϺ����ս���ṹ��
fprintf('�˲�...%5.0f %%',0);
[M,N]=size(SINS);%Ԥ����������SINS�ӿ� M=1Ϊ��SINS��
gvar;    % ����ȫ�ֱ���
looplen=floor(len);
ts=SINS{1,1}.ts;%����ʱ��
nts=2*ts;%����ʱ��
%��ʼ����
Rloop=cell(M,N);%��Ե��� R ���Ʊ���
Uloop=cell(M,N);%��Ե��� U ���Ʊ���
INSSFR=cell(M,N);%��Ե����˲����
Filter=cell(M,N);%��Ե����˲��м�ֵ
vn=cell(M,N);%��Ե����˲��ٶ��м�ֵ
qms=cell(M,N);%��Ե��� q ���Ʊ���
attnow=cell(M,N);%�˲� �Ƕ� ���Ʊ���
dR=cell(M,N);%�˲����� dR ���Ʊ���
attint=cell(M,N);%��ʼ�� ��������Ƕ���� ��ת����yxz
kfft=cell(M,N);%�˲�������״̬ת�ƾ���
kf=cell(M,N);%�������˲���
RFilter=cell(M,N);%��ͨ�˲�����ʼֵ
ws_m=cell(M,N);%ѭ��ʱÿ���ӹߵ��������� ���ٶ�
fs_m=cell(M,N);%ѭ��ʱÿ���ӹߵ��������� ����

NINSS.uf=zeros(3,looplen);
for m=1:M
    for n=1:N
        Rloop{m,n}=SINS{m,n}.R0;
        Uloop{m,n}=cross(MINS.wim0,SINS{m,n}.R0);

        INSSFR{m,n}.Rall=zeros(looplen,3);
        INSSFR{m,n}.vnall=zeros(looplen,3);
        INSSFR{m,n}.qmsall=zeros(looplen,4);
        INSSFR{m,n}.attall=zeros(3,looplen);

        Filter{m,n}.X=zeros(15,looplen);
        
        Cmserr_1=my_a2mat(atterr{m,n});%ע���м䴦����ת˳��yxz ��׼ȷ�ĳ�ʼ��
        qms{m,n}=m2qua(Cmserr_1);%��ʼ��̬ ��׼ȷ�ĳ�ʼ��
        attint{m,n}=m2att(Cmserr_1);%��ת����yxz ������θ����� ����� ��λ��

        fs0=[0;0;0];
        kfft{m,n}=my_kfft15(MINS.wim0,q2mat(qms{m,n})',fs0,nts);
        kf{m,n}= kfinit(KFinit{m,n}.Qk, KFinit{m,n}.Rk, KFinit{m,n}.P0,kfft{m,n}.phi,kfft{m,n}.H);  % kf�˲�����ʼ��
        kf{m,n}.Xk(10:12) = [0.1;0.1;0.1]*dph;
        kf{m,n}.Xk(13:15) = [20;20;20]*ug;
        RFilter{m,n}=SINS{m,n}.R0;
    end
end

for i=1:looplen
    if flag.select==0%�����
        for m=1:M
            for n=1:N
                ws_m{m,n}=SINS{m,n}.wis(:,(2*i-1):(2*i))';
                fs_m{m,n}=SINS{m,n}.fs(:,(2*i-1):(2*i))';
            end
        end
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; 
        fm_m=MINS.fm(:,(2*i-1):(2*i))';
    elseif flag.select==1%ֱ��ʹ����Ե��������
        SINSR=selcCell;
        for m=1:M
            for n=1:N
            ws_m{m,n}=SINSR{m,n}.ws_m_addnoise(:,:,i);
            fs_m{m,n}=SINSR{m,n}.fs_m_addnoise(:,:,i);
            end
        end
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; 
        fm_m=MINS.fm(:,(2*i-1):(2*i))';
    else %��ע�����
        SERR=selcCell;
        for m=1:M
            for n=1:N
                [ws_m{m,n}, fs_m{m,n}] = my_imuadderr...
                    (SINS{m,n}.wis(:,(2*i-1):(2*i))',...
                    SINS{m,n}.fs(:,(2*i-1):(2*i))', ...
                    SERR{m,n}.eb,SERR{m,n}.web,...
                    SERR{m,n}.db,SERR{m,n}.wdb,ts);%�ӹ�ע������
            end
        end
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; 
        fm_m=MINS.fm(:,(2*i-1):(2*i))';
    end
%��������������������������������������������������������������������������������������������������
% ���¿��Ի��ÿ��ʱ�̵����ߣ�1�������ӹߣ�M*N��������������
% ����Ϊ wm_m �������߽����� fm_m �������߱���
%        ws_m �����ӹ߽����� fs_m �����ӹ߱��� ��Ԫ��{M,N}ʹ��
% ����Ϊ�˺�����ڵ�����������ʹ��
    for m=1:M
        for n=1:N
            [qms{m,n},vn{m,n},Uloop{m,n},Rloop{m,n}] = my_relinsupdate5...
                (qms{m,n},Uloop{m,n},Rloop{m,n},wm_m,fm_m,ws_m{m,n},fs_m{m,n},ts);% �ӹߴ���Ե���
            RFilter{m,n}=my_Ofilter(RFilter{m,n},Rloop{m,n},2,nts);%��ͨ�˲�
        end
    end
    if flag.EnFusion  
        for m=1:M
            for n=1:N
                attnow{m,n}=q2att(qms{m,n});
                uf=attnow{m,n}-attint{m,n};%����

                dR{m,n}=my_getdR(RFilter{m,n},SINS{m,n}.R0,uf);%�������

                wim=sum(wm_m)/2;%����������
                fs=sum(fm_m)/2;
                kfft{m,n}=my_kfft15(wim',q2mat(qms{m,n})',fs',nts);
                kf{m,n}.Phikk_1=kfft{m,n}.phi;
                kf{m,n}.Gammak=kfft{m,n}.Gammak;

                kf{m,n} = kfupdate(kf{m,n},dR{m,n},'B');%�������˲�

                Filter{m,n}.X(:,i)=kf{m,n}.Xk;
                
                attint{m,n}=q2att(qdelphi(a2qua(attint{m,n}),-kf{m,n}.Xk(1:3)));% 2021 03 15���� ���³�ʼ��̬�����
                
                qms{m,n} = qdelphi(qms{m,n},-kf{m,n}.Xk(1:3));
                kf{m,n}.Xk(1:3) = 0;  % ������
                Uloop{m,n}= Uloop{m,n}- kf{m,n}.Xk(4:6);  kf{m,n}.Xk(4:6) = 0;
                Rloop{m,n}= Rloop{m,n}- kf{m,n}.Xk(7:9);  kf{m,n}.Xk(7:9) = 0;

                INSSFR{m,n}.qmsall(i,:)=qms{m,n}';
                INSSFR{m,n}.Rall(i,:)=Rloop{m,n}';
                INSSFR{m,n}.vnall(i,:)=vn{m,n}';

                [rz,rx,ry]=dcm2angle(q2mat(INSSFR{m,n}.qmsall(i,:)),'zxy');%ת������ͳ������
                INSSFR{m,n}.attall(:,i)=[rx,ry,-rz];
            end
        end
        qms_all=my_CellAvr(qms);
        Uloop_all=my_CellAvr(Uloop);
        vn_all=my_CellAvr(vn);
        Rloop_all=my_CellAvr(Rloop);
        Rfall=Rloop_all-SINS{m,n}.R0;
        if Rfall(3)>0
            NINSS.uf(:,i)=[0;asin(sqrt(abs(Rfall'*Rfall)/abs(SINS{m,n}.R0'*SINS{m,n}.R0)));0];
        else
            NINSS.uf(:,i)=[0;-asin(sqrt(abs(Rfall'*Rfall)/abs(SINS{m,n}.R0'*SINS{m,n}.R0)));0];
        end
        
        for m=1:M % ��С���˷��� û�и��·�����
            for n=1:N
%                 kf{m,n}.Xk(1:3) = q2rv(qdelphi(qms{m,n},qms_all));
%                 qms{m,n}= qms_all;  
%                 kf{m,n}.Xk(4:6) = Uloop{m,n}-Uloop_all;
%                 Uloop{m,n}= Uloop_all;
                kf{m,n}.Xk(7:9) = Rloop{m,n}-Rloop_all;
                Rloop{m,n}= Rloop_all;
            end
        end
        
        NINSS.qmsall(i,:)=qms_all';
        NINSS.Rall(i,:)=Rloop_all';
        NINSS.vnall(i,:)=vn_all';
        attMat=a2mat(2*NINSS.uf(:,i));
        [rz,rx,ry]=dcm2angle(attMat,'zxy');%ת������ͳ������
        NINSS.attnewall(:,i)=[rx,ry,-rz];
        [rz,rx,ry]=dcm2angle(q2mat(NINSS.qmsall(i,:)),'zxy');%ת������ͳ������
        NINSS.attall(:,i)=[rx,ry,-rz];
        
    else
        
        for m=1:M
            for n=1:N
                attnow{m,n}=q2att(qms{m,n});
                uf=attnow{m,n}-attint{m,n};%����

                dR{m,n}=my_getdR(RFilter{m,n},SINS{m,n}.R0,uf);%�������

                wim=sum(wm_m)/2;%����������
                fs=sum(fm_m)/2;
                kfft{m,n}=my_kfft15(wim',q2mat(qms{m,n})',fs',nts);
                kf{m,n}.Phikk_1=kfft{m,n}.phi;
                kf{m,n}.Gammak=kfft{m,n}.Gammak;

                kf{m,n} = kfupdate(kf{m,n},dR{m,n},'B');%�������˲�

                Filter{m,n}.X(:,i)=kf{m,n}.Xk;

            %     [rz,rx,ry]=dcm2angle(q2mat(rv2q(kf{m,n}.Xk(1:3))),'zxy');
            %     Filter{m,n}.atterr=[rx;ry;-rz];

            %     kf{m,n}.Xk(1:3) = 0;
                qms{m,n} = qdelphi(qms{m,n},-kf{m,n}.Xk(1:3));
                kf{m,n}.Xk(1:3) = 0;  % ������
                Uloop{m,n}= Uloop{m,n}- kf{m,n}.Xk(4:6);  kf{m,n}.Xk(4:6) = 0;
                Rloop{m,n}= Rloop{m,n}- kf{m,n}.Xk(7:9);  kf{m,n}.Xk(7:9) = 0;

                INSSFR{m,n}.qmsall(i,:)=qms{m,n}';
                INSSFR{m,n}.Rall(i,:)=Rloop{m,n}';
                INSSFR{m,n}.vnall(i,:)=vn{m,n}';

                [rz,rx,ry]=dcm2angle(q2mat(INSSFR{m,n}.qmsall(i,:)),'zxy');%ת������ͳ������
                INSSFR{m,n}.attall(:,i)=[rx,ry,-rz];
            end
        end
        
        NINSS=INSSFR{1,1};%�����ںϷ���ֵ ��ʱ���ں�
    end
%����������������������������������������������������������������������������������������������������������������
    if mod(i,looplen/100)==0
              fprintf('\b\b\b\b\b\b\b%5.0f %%', i/looplen*100);			%������ʾ
    end
end
fprintf('\n�˲���ɣ�\n'); 