function err=my_getSErr(ex,len,dir)
% �������ƣ�my_getSErr
% �������ܣ��������յ����
% ���룺ex      :��ֵ
%       dir     :����1�����м��㣬2���м���
% �����INSSFR:��Ե������
len=floor(len);
if dir==1
    err=sum(ex(end-len+1:end,:),dir)/len;
elseif dir==2
    err=sum(ex(:,end-len+1:end),dir)/len;
end