function fe=my_fe1(I0,n,sita)
% �벨�������߹�һ������ͼ����
% fe=I0*n/(2*pi)*(cos(pi/2.*cos(sita)))./sin(sita);
% �벨���ӵĵ�Ԫ����ͼ����
fe=I0*n/(2*pi)*(sin(pi/2.*cos(sita+pi/2)))./cos(sita+pi/2);