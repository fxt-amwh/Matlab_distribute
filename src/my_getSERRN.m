function SERR=my_getSERRN(SERRList)
% �������ƣ�KFinitN
% �������ܣ����󿨶����˲���
% ���룺SERR      :�ӹ����Ԫ��
%       SINS      :�ӹ���ϢԪ��  
%       nts = 2*ts;�������Ͳ���ʱ��
% �����atterr:��Ե������
[~,N,M]=size(SERRList.eb);%Ԥ����������SINS�ӿ� M=1Ϊ��SINS��
SERR=cell(M,N);
for m=1:M
    for n=1:N
        s.eb = SERRList.eb(:,n,m); s.web = SERRList.web(:,n,m);   %���ݳ�ֵ��ƫ���Ƕ��������
        s.db = SERRList.db(:,n,m); s.wdb = SERRList.wdb(:,n,m);  %���ٶȼƳ�ֵƫֵ���ٶ��������
        SERR{m,n}=s;
    end
end
