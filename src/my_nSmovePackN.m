function SmoveCell=my_nSmovePackN(Sinf)
% �������ƣ�my_nSmovePackN
% �������ܣ���������Ե����ӹ����������ݽ�� ����
% ���룺Sinf      :�ӹ���Ϣ���þ��� Ԫ��
% �����SmoveCell :�ӹ߽ṹ��Ԫ��
M=size(Sinf,2);%Ԥ����������SINS�ӿ� M=1Ϊ��SINS��
N=size(Sinf{1}.Rlist,2);
SmoveCell=cell(M,N);
for m=1:M
    SmoveCellAllN=my_nSmovePack(Sinf{m});
    for n=1:N
        SmoveCell{m,n}=SmoveCellAllN{1,n};
    end
end