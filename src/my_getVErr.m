function err=my_getVErr(ex,beg,End,dir)
% �������ƣ�my_getVErr
% �������ܣ��������յľ�����
% ���룺beg      :��ֵ
%       end     :����1�����м��㣬2���м���
%       dir     :����1�����м��㣬2���м���
% �����INSSFR:��Ե������
if dir==1
    err=sqrt(var(ex(beg:End,:),0,dir));
elseif dir==2
    err=sqrt(var(ex(:,beg:End),0,dir));
end