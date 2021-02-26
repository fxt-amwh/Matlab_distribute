function [qms,vn,R,qnb_m, vn_m, pos_m ,qnb_s, vn_s, pos_s] = my_relinsupdate3(qms,vn,R,qnb_m, vn_m, pos_m, wm_m, vm_m,qnb_s, vn_s, pos_s, wm_s, vm_s, ts)
nn = size(wm_m,1);  nts = nn*ts;
%% ����
[phim, dvbm] = cnscl(wm_m, vm_m);  % Բ׶���/��������
eth = earth(pos_m, vn_m);  % ������ز�������
vn1_m = vn_m + rv2m(-eth.wnin*nts/2)*qmulv(qnb_m,dvbm) + eth.gcc*nts;  % �ٶȸ���
vn_ma = (vn_m+vn1_m)/2;
pos_m = pos_m + [vn_ma(2)/eth.RMh;vn_ma(1)/eth.clRNh;vn_ma(3)]*nts;  vn_m = vn1_m;  % λ�ø���
qnb_m = qmul(rv2q(-eth.wnin*nts), qmul(qnb_m, rv2q(phim)));  % ��̬����
qnb_m = qnormlz(qnb_m);
%% �ӹ�
[phim, dvbm] = cnscl(wm_s, vm_s);  % Բ׶���/��������
eths = earth(pos_s, vn_s);  % ������ز�������
vn1_s = vn_s + rv2m(-eths.wnin*nts/2)*qmulv(qnb_s,dvbm) + eths.gcc*nts;  % �ٶȸ���
vn_sa = (vn_s+vn1_s)/2;

pos_s = pos_s + [vn_sa(2)/eths.RMh;vn_sa(1)/eths.clRNh;vn_sa(3)]*nts;  vn_s = vn1_m;  % λ�ø���
qnb_s = qmul(rv2q(-eths.wnin*nts), qmul(qnb_s, rv2q(phim)));  % ��̬����
qnb_s = qnormlz(qnb_s);

vn_msa=(q2mat(qnb_s)')*vn_sa-(q2mat(qnb_m)')*vn_ma;

R=R+vn_msa*nts;
phim = my_getdphim(wm_m);
phis = my_getdphim(wm_s);
qms = qmul(rv2q(-phis), qmul(qms, rv2q(phim)));  % ��̬����
qms = qnormlz(qms);





