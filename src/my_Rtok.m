function k=my_Rtok(R,R0)
% �������ƣ�my_getdR
% �������ܣ���R���������R�ĵ�λ����
% ���룺R      :��ǰ�ӹ�����������ϵλ��
%     :R0      :��ǰ�ӹ߱��λ��
% �����k:λ���������ֵ
uf=my_Rtouf(R,R0);
k=[0;uf(2)/sqrt(R0'*R0);0];
% my_getdRf(L,u,flag)


