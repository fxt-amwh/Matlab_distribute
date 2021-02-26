function att=my_q2attzxy(q)
[rz,rx,ry]=dcm2angle(q2mat(q),'zxy');%转换到传统的坐标
att=[rx,ry,-rz]';
