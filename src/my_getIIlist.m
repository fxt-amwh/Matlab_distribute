function anti=my_getIIlist(dim,T)
% �������ƣ�my_getIIlist
% �������ܣ�ͨ��ԭ�����ڲ��������ڵ�ֵ��ñ�������ֵ
% ���룺dim:��ǰԭ�����ڵ��β���ʱ��������
%      T  : ����ʱ��
% �����anti:���������ڵ�ԭ����ƽ��ֵ����
n=length(dim);
anti=zeros(1,n);
for i=2:n
    anti(i-1)=(dim(i)-dim(i-1))/T;
end