function uf=my_Rtouf(R,R0)
% �������ƣ�my_getdR
% �������ܣ���R���������
% ���룺R      :��ǰ�ӹ�����������ϵλ��
%     :R0      :��ǰ�ӹ߱��λ��
% �����uf:λ���������ֵ
Rf=R-R0;
if Rf(3)>0
    uf=[0;asin(sqrt((Rf'*Rf)/(R0'*R0)));0];
else
    uf=[0;-asin(sqrt((Rf'*Rf)/(R0'*R0)));0];
end
