function [qnb, vn, pos, eth] = my_SINSVgetResult(qnb, vn, pos, wm, vm, ts)
% �������ƣ�my_SINSVgetResult
% �������ܣ����ӹ����������Ե���
% ���룺MINS     �����߽ṹ��
%       SINS     ���ӹ߽ṹ��
%       varargin:������΢�����������,��Ԫ�����          
%              ÿһ�����������δ��Cms��������������任����Csm���֣�3*3
%                                  dU  �ٶ�΢��
%                                  fm  ���߱���
%                                  wim ��������ڹ�������ϵ�ٶ�����������ϵ��ͶӰ
% �����INSSR:��Ե������
nn = size(wm,1);  nts = nn*ts;
[phim, dvbm] = cnscl(wm, vm);  % Բ׶���/��������
eth = earth(pos, vn);  % ������ز�������
vn1 = vn + rv2m(-eth.wnin*nts/2)*qmulv(qnb,dvbm);  % �ٶȸ���
vn = (vn+vn1)/2;
pos = pos + [vn(2)/eth.RMh;vn(1)/eth.clRNh;vn(3)]*nts;  vn = vn1;  % λ�ø���
qnb = qmul(rv2q(-eth.wnin*nts), qmul(qnb, rv2q(phim)));  % ��̬����
qnb = qnormlz(qnb);