function INSSR=my_SINSgetResultN(MINS,SINS,len,SERR,atterr)
% �������ƣ�my_SINSgetResultN
% �������ܣ����ӹ����������Ե���
% ���룺MINS   :���ߵ�
%       SINS   :�ӹߵ�����Ԫ�� M*N 
%       len    :�������г��ȵ�һ�� 1*1
%       SERR   :�ӹߵ����Ԫ�� M*N 
%       atterr :�ӹߵ�λ�����Ԫ�� M*N 
% �����INSSR:��Ե������ ��ӿ�
%% ��Ե���
fprintf('���ӹ���Ե���...%5.0f %%',0); 
[M,N]=size(SINS);%Ԥ����������SINS�ӿ� M=1Ϊ��SINS��
INSSR=cell(M,N);
for m=1:M
    for n=1:N
        INSSR{m,n}=my_SINSgetResult(MINS,SINS{m,n},len,SERR{m,n},atterr{m,n});
    end
end
fprintf('\n������Ե�����ɣ�\n'); 