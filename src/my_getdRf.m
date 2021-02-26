function Rf=my_getdRf(L,u,flag)
% �������ƣ�my_getdRf
% �������ܣ���Ե����˲�״̬ת�ƾ���
% ���룺L      :���λ��     
%       u      :������ν�
% �����Rf:�ӹ�λ��
if nargin<3
    flag=1;
end
ud2=u/2;
switch flag
    case 1%С�ǶȽ���
        Rf=[-L(1)*(ud2(2)^2)+L(2)*ud2(3)
            -L(2)*(ud2(3)^2)+L(3)*ud2(1)
            -L(3)*(ud2(1)^2)+L(1)*ud2(2)];
    case 2%��ȷ����
        Rf=[-L(1)*sin(ud2(2))*sin(ud2(2))+L(2)*sin(ud2(3))*cos(ud2(3))
            -L(2)*sin(ud2(3))*sin(ud2(3))+L(3)*sin(ud2(1))*cos(ud2(1))
            -L(3)*sin(ud2(1))*sin(ud2(1))+L(1)*sin(ud2(2))*cos(ud2(2))];
end