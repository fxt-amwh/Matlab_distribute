lenmax=2;
xdir=linspace(0,lenmax,1000);
[~,dz]=my_fixdefload(xdir,lenmax,10*pi/180,1);%�����غ�����ģ�ͱ�׼�� ��׼ģ��1
[~,qD]=my_xtodz(0,lenmax,15*pi/180);%�������λ�á���������ȡ�غ�ϵ��
qD
figure
plot(xdir,dz)
ylim([0,lenmax])
xstd=0.5;
[ddz,ddx,i]=my_XstdToDz(xstd,2,5*pi/180)%���5������ʱ���
[ddz,ddx,i,a]=my_XstdToDzAndAtt(xstd,2,5*pi/180)%���5������ʱ���
a*180/pi
atan(ddz/(xstd+ddx))*180/pi
[ddz,ddx,i,a,astd]=my_XstdToDzAndAtt(xstd,2,5*pi/180)%���5������ʱ���
a*180/pi
astd*180/pi