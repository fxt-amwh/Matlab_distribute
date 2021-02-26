function [phim] = my_getdphim(wm)
cs = [  [2,    0,    0,    0,   0    ]/3
      [9,    27,   0,    0,    0    ]/20
      [54,   92,   214,  0,    0    ]/105
      [250,  525,  650,  1375, 0    ]/504
      [2315, 4558, 7296,  7834, 15797]/4620  ];  % 2-6子样补偿系数
wmm = sum(wm,1);
n = size(wm, 1);  % 子样数
if n>1
    csw = cs(n-1,1:n-1)*wm(1:n-1,:); 
    dphim = cross(csw,wm(n,:));  % 圆锥补偿量
end
phim = (wmm+dphim)';