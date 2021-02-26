%% 信号低通滤波器测试
clc
clear
close all
ts=0.001;
t=0:ts:600;
len=length(t);
f1=10;%高频分量
f2=0.01;%低频分量
y1=5*sin(2*pi*f1*t);
y2=100*sin(2*pi*f2*t);
y=y1+y2;
% fc=10;
% Tc=1/(2*pi*fc);
% a=ts/(Tc+ts);
ONAN=zeros(1,len);
O10=zeros(1,len);
O2=zeros(1,len);
O002=zeros(1,len);
ONAN(1)=my_Ofilter(0,y(1),1);%a=1表示带宽频率无穷
O10(1)=my_Ofilter(0,y(1),10,ts);%10Hz
O2(1)=my_Ofilter(0,y(1),2,ts);%2Hz
O002(1)=my_Ofilter(0,y(1),0.02,ts);%0.02Hz
for i=2:len
    O10(i)=my_Ofilter(O10(i-1),y(i),10,ts);
    ONAN(i)=my_Ofilter(ONAN(i-1),y(i),1);%a=1表示带宽频率无穷
    O10(i)=my_Ofilter(O10(i-1),y(i),10,ts);%10Hz
    O2(i)=my_Ofilter(O2(i-1),y(i),2,ts);%2Hz
    O002(i)=my_Ofilter(O002(i-1),y(i),0.02,ts);%0.02Hz
end
figure
subplot(221)
plot(t,y,t,ONAN);
legend('原始信号','NAN滤波结果');
subplot(222)
plot(t,y,t,O10);
legend('原始信号','10Hz滤波结果');
subplot(223)
plot(t,y,t,O2);
legend('原始信号','2Hz滤波结果');
subplot(224)
plot(t,y,t,O002);
legend('原始信号','0.02HZ滤波结果');