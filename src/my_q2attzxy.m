function att=my_q2attzxy(q)
[rz,rx,ry]=dcm2angle(q2mat(q),'zxy');%ת������ͳ������
att=[rx,ry,-rz]';
