function uf=my_Rftouf(Rf,R0)
% �������ƣ�my_getdR
% �������ܣ���R���������
% ���룺Rf      :��ǰ�ӹ�����������ϵλ��
%     :R0      :��ǰ�ӹ߱��λ��
% �����uf:λ���������ֵ
if Rf(3)>0
    uf=[0;asin(sqrt((Rf'*Rf)/(R0'*R0)));0];
else
    uf=[0;-asin(sqrt((Rf'*Rf)/(R0'*R0)));0];
end