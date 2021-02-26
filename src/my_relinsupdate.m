function [qms, vn, R] = my_relinsupdate(qms, vn, R, wm_m, vm_m,ws_m,vs_m,ts)
nts = size(wm_m,1)*ts;%���θ���ʱ����
Csm=q2mat(qms)';
dv=my_relgetdv(Csm,vn',R',wm_m, vm_m,ws_m,vs_m,ts)';
vn1 = vn + dv;  % �ٶȸ���
vn = (vn+vn1)/2;
R = R + vn*nts;  vn = vn1;  % λ�ø���
phim = my_getdphim(wm_m);
phis = my_getdphim(ws_m);
qms = qmul(rv2q(-phis), qmul(qms, rv2q(phim)));  % ��̬����
qms = qnormlz(qms);

