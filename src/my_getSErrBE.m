function err=my_getSErrBE(ex,beg,End,dir)
% �������ƣ�my_getSErr
% �������ܣ��������յ����
% ���룺ex      :��ֵ
%       dir     :����1�����м��㣬2���м���
% �����INSSFR:��Ե������
len=End-beg+1;
if dir==1
    err=sum(abs(ex(beg:End,:)),dir)/len;
elseif dir==2
    err=sum(abs(ex(:,beg:End)),dir)/len;
end