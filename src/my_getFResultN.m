function [INSSFR,Filter]=my_getFResultN(MINS,SINS,KFinit,atterr,len,selcCell,select)
% �������ƣ�my_SINSgetFResultN
% �������ܣ����ӹ����������Ե����˲�
% ���룺MINS      :����
%       SINS      :�ӹ�
%       KFinit    :�ӹ��˲�����ʼֵ
%       atterr    :�ӹ߰�װ���ǹ��Ƴ�ʼֵ
%       len       :�˲�����       
%       selcCell  :select=1��selcCell��ʾ SINSR �ӹߴ���Ե������
%                  select=2��selcCell��ʾ SERR  �ӹ�����ʼֵ
%       select    :ѡ���ӹ������Դ
% �����INSSFR:��Ե������
fprintf('�˲�...%5.0f %%',0);
[M,N]=size(SINS);%Ԥ����������SINS�ӿ� M=1Ϊ��SINS��
INSSFR=cell(M,N);
Filter=cell(M,N);
if select==1%ֱ��ʹ����Ե��������
    SINSR=selcCell;
    for m=1:M
        for n=1:N
        [INSSFR{m,n},Filter{m,n}]=my_getFResult...
            (MINS,SINS{m,n},KFinit{m,n},atterr{m,n},len,...
        SINSR{m,n}.ws_m_addnoise,SINSR{m,n}.fs_m_addnoise);
        end
    end
else%����ע�����
    SERR=selcCell;
    for m=1:M
        for n=1:N
        [INSSFR{m,n},Filter{m,n}]=my_getFResult...
            (MINS,SINS{m,n},KFinit{m,n},atterr{m,n},len,SERR{m,n});
        end
    end
end
fprintf('\n�˲���ɣ�\n'); 

