function atterr=my_getAtterr(aList)
% �������ƣ�KFinitN
% �������ܣ����󿨶����˲���
% ���룺SERR      :�ӹ����Ԫ��
%       SINS      :�ӹ���ϢԪ��  
%       nts = 2*ts;�������Ͳ���ʱ��
% �����atterr:��Ե������
[~,N,M]=size(aList);%Ԥ����������SINS�ӿ� M=1Ϊ��SINS��
atterr=cell(M,N);
for m=1:M
    for n=1:N
        att=aList(:,n,m);
        atterr{m,n}=att;
    end
end
