clc
close all
g2d=180/pi;
d2g=1/g2d;
uf_lw_m=5*d2g;
f=10;
ts=0.01;
T=100;
t=0:ts:(T/f);
len=length(t);
uf_lw=[zeros(1,len);uf_lw_m*sin(2*pi*t);zeros(1,len)];
R0=[15;0;0];
Rf0=[0;0;0];%理想时 子0时刻前偏移
Rf=[-R0(1)*sin(uf_lw(2,:)).*sin(uf_lw(2,:));zeros(1,len);R0(1)*sin(uf_lw(2,:)).*cos(uf_lw(2,:))];
R=R0+Rf;
figure
subplot(311)
plot(t,R(1,:));
subplot(312)
plot(t,R(2,:));
subplot(313)
plot(t,R(3,:));

l_list=[linspace(0,R0(1),5);zeros(2,5)];
N=size(l_list,2);
Rf0=zeros(3,N);
Rf=zeros(3,len,N);
R=zeros(3,len,N);
now_t=1;
for i=1:N
    k=l_list(1,i)/R0(1);
    Rf(:,:,i)=[-l_list(1,i)*sin(k*uf_lw(2,:)).*sin(k*uf_lw(2,:));zeros(1,len);l_list(1,i)*sin(k*uf_lw(2,:)).*cos(k*uf_lw(2,:))];
    R(:,:,i)=l_list(:,i)+Rf(:,:,i);
end
figure
subplot(311)
plot(t,R(1,:,3));
subplot(312)
plot(t,R(2,:,3));
subplot(313)
plot(t,R(3,:,3));
figure
plot(t,R(3,:,3));
figure
for i=1:len
    x1(:,:)=R(1,now_t,:);
    p1(:,:)=R(3,now_t,:);
    plot(x1,p1,'-*');
    axis([0,15,-15/2,15/2])
    pause(ts);
    now_t=now_t+1;
end



    