function anti=my_getII(dim,T)
% �������ƣ�my_getII
% �������ܣ�ͨ�������ڲ��������ڵĻ��ֻ�ñ�������ֵ
% ���룺dim:��ǰ���������ڵ��β���ʱ���ڵĻ�������ֵ
%      T  : ����ʱ��
% �����anti:���������ڵı�������ƽ��ֵ����
if length(dim)==length(T)
    anti=dim./T;
else
    anti=dim/T(1);
end
