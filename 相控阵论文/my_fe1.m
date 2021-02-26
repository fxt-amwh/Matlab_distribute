function fe=my_fe1(I0,n,sita)
% 半波振子天线归一化方向图函数
% fe=I0*n/(2*pi)*(cos(pi/2.*cos(sita)))./sin(sita);
% 半波振子的单元方向图函数
fe=I0*n/(2*pi)*(sin(pi/2.*cos(sita+pi/2)))./cos(sita+pi/2);