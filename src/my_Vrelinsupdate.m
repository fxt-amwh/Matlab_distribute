function [qnm,vnm,posm,qns,vns,poss,VL,ethm] = my_Vrelinsupdate(qnm,vnm,posm,wm,vm,qns,vns,poss,ws,vs,R,ts)
% �������ƣ�my_Vrelinsupdate
% �������ܣ����ӹ����������Ե���
% ���룺qnm     ��������̬��Ԫ�� ע�⺬�壬������ʷԭ�������qnm��ʾ����������ϵ����������ϵ��Ԫ��
%       vnm     �������ٶ���� 
%       posm    �����߾�γ��
%       wm      ��������������
%       vm      �������ٶ�����
%       qns     ���ӹ���̬��Ԫ�� ע�⺬�壬������ʷԭ�������qns��ʾ���ӹ�����ϵ����������ϵ��Ԫ��
%       vns     ���ӹ��ٶ����
%       poss    ���ӹ߾�γ��
%       ws      ���ӹ���������
%       vs      ���ӹ��ٶ�����
%       R       �����Ӽ�˱�
%       ts      ������ʱ��
% �����qnm     ��������̬��Ԫ�� ע�⺬�壬������ʷԭ�������qnm��ʾ����������ϵ����������ϵ��Ԫ��
%       vnm     �������ٶ����
%       posm    �����߾�γ��
%       qns     ���ӹ���̬��Ԫ�� ע�⺬�壬������ʷԭ�������qns��ʾ���ӹ�����ϵ����������ϵ��Ԫ��
%       vns     ���ӹ��ٶ����
%       poss    ���ӹ߾�γ��
%       VL      ����������ϵ���ӹ߸˱��ٶ�
%       ethm    �����߱�׼�µ�����ز���
%����
nn = size(wm,1);  nts = nn*ts;
[phim, dvbm] = cnscl(wm, vm);  % Բ׶���/��������
ethm = earth(posm, vnm);  % ������ز�������
vn1 = vnm + rv2m(-ethm.wnin*nts/2)*qmulv(qnm,dvbm);  % �ٶȸ���
vnm = (vnm+vn1)/2;
posm = posm + [vnm(2)/ethm.RMh;vnm(1)/ethm.clRNh;vnm(3)]*nts;  vnm = vn1;  % λ�ø���
qnm = qmul(rv2q(-ethm.wnin*nts), qmul(qnm, rv2q(phim)));  % ��̬����
qnm = qnormlz(qnm);
%�ӹ�
[phis, dvbs] = cnscl(ws, vs);  % Բ׶���/��������
eths = earth(poss, vns);  % ������ز�������
vns1 = vns + rv2m(-eths.wnin*nts/2)*qmulv(qns,dvbs);  % �ٶȸ���
vns = (vns+vns1)/2;
poss = poss + [vns(2)/eths.RMh;vns(1)/eths.clRNh;vns(3)]*nts;  vns = vns1;  % λ�ø���
qns = qmul(rv2q(-eths.wnin*nts), qmul(qns, rv2q(phis)));  % ��̬����
qns = qnormlz(qns);
%�˱��ٶ�
w_MiM=mean(wm,1)/ts;%������ת��Ϊ���ٶ�
w_nie=ethm.wnie;
CnM=q2mat(qnm);
% VL=CMn*cross(w_MiM,R)'-cross(w_nie,CMn*R);
% VL=CMn*(askew(w_MiM'-CnM*w_nie))*R;
VL=CnM*cross(w_MiM,R)'-askew(w_nie)*CnM*R;