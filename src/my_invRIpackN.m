function INSS=my_invRIpackN(MINS,Smove)
% �������ƣ�my_invRIpackN
% �������ܣ���������Ե����ӹ����������ݽ�� ���Ӳ�ȷ������
% ���룺MINS      :����
%       Smove     :�ӹ��˶�Ԫ�� Ԫ��������Ӧ�ӹ�����
% �����INSS      :�ӹ߽ṹ��Ԫ��
[M,N]=size(Smove);%Ԥ����������SINS�ӿ� M=1Ϊ��SINS��
INSS=cell(M,N);
for m=1:M
    for n=1:N
        INSS{m,n}=my_invRIpack(MINS,Smove{m,n});
    end
end
