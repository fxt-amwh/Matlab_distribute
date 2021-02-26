function [wm] = att2wm(att,att0)%惯性传感器信息生成（反演算法）
wm0 = zeros(3,1); I33 = eye(3);
wm = att(:,1:end);
qbb = qmul(qconj(a2qua(att0)),a2qua(att(:,1)));
phim = q2rv(qbb);
wm1 = (I33+askew(1/12*wm0))\phim;
wm(:,1) = wm1;   wm0 = wm1;
for k=2:length(att)
 
    qbb = qmul(qconj(a2qua(att(:,k-1))),a2qua(att(:,k)));
    phim = q2rv(qbb);
    wm1 = (I33+askew(1/12*wm0))\phim;
    wm(:,k-1) = wm1;   wm0 = wm1;
end