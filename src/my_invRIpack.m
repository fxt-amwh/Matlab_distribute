function INSS=my_invRIpack(MINS,Smove)
% �������ƣ�my_invRIpack
% �������ܣ���������Ե����ӹ����������ݽ��
% ���룺MINS      :����
%       Smove     :�ӹ��˶�
% �����INSS      :�ӹ߽ṹ��
INSS=my_invRI(Smove.R,Smove.att,Smove.atterr0,...
    MINS.wim,MINS.fm,MINS.ts,Smove.R0,Smove.att_s0,Smove.U0);
