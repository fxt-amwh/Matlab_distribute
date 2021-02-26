lenmax=2;
xdir=linspace(0,lenmax,1000);
[~,dz]=my_fixdefload(xdir,lenmax,10*pi/180,1);%机翼载荷挠曲模型标准量 精准模型1
[~,qD]=my_xtodz(0,lenmax,15*pi/180);%根据翼尖位置、挠曲角求取载荷系数
qD
figure
plot(xdir,dz)
ylim([0,lenmax])
xstd=0.5;
[ddz,ddx,i]=my_XstdToDz(xstd,2,5*pi/180)%获得5°挠曲时误差
[ddz,ddx,i,a]=my_XstdToDzAndAtt(xstd,2,5*pi/180)%获得5°挠曲时误差
a*180/pi
atan(ddz/(xstd+ddx))*180/pi
[ddz,ddx,i,a,astd]=my_XstdToDzAndAtt(xstd,2,5*pi/180)%获得5°挠曲时误差
a*180/pi
astd*180/pi