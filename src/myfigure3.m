if flag.Ceil==0
    R=SINS1.R;%理论位置
    RFR=SINSFR1.Rall(1:2:end,:);
    T=ts:ts:(Sinf.len*ts);
else
    R=SINS{flag.flagSINS_M,flag.flagSINS_N}.R;%理论位置
    RFR=SINSFR{flag.flagSINS_M,flag.flagSINS_N}.Rall(1:2:end,:);
    T=ts:ts:(SinfCell{1,1}.len*ts);
end
figure%理论相对位置与滤波解算相对位置对比
subplot(311)
plot(T(2:4:end),R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,1),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('x/m');
title('主子惯导相对位置R（滤波）');
grid on;
subplot(312)
plot(T(2:4:end),R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,2),'--','LineWidth',2)
grid on;
legend('真实值','计算值');
xlabel('时间/s');ylabel('y/m');
title('主子惯导相对位置R（滤波）');
subplot(313)
plot(T(2:4:end),R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,3),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('z/m');
title('主子惯导相对位置R（滤波）');
grid on;