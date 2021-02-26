if flag.Ceil==0
    R=SINS1.R;%理论位置
    RFR=NSINSFR.Rall(1:2:end,:);
    T=ts:ts:(Sinf.len*ts);
else
    R=SINS{flag.flagSINS_M,flag.flagSINS_N}.R;%理论位置
    RFR=NSINSFR.Rall(1:2:end,:);
    T=ts:ts:(SinfCell{1,1}.len*ts);
end
figure
subplot(311)
plot(T(2:4:end),RFR(:,1)'-R(1,2:4:end),'LineWidth',2)
legend('X');
xlabel('时间/s');ylabel('x/m');
title('天线中心相对位置误差');
grid on;ylim([-0.01,0.01]);
subplot(312)
plot(T(2:4:end),RFR(:,2)'-R(2,2:4:end),'LineWidth',2)
grid on;
legend('Y');
xlabel('时间/s');ylabel('y/m');
title('天线中心相对位置误差');ylim([-0.01,0.01]);
subplot(313)
plot(T(2:4:end),RFR(:,3)'-R(3,2:4:end),'LineWidth',2)
legend('Z');
xlabel('时间/s');ylabel('z/m');
title('天线中心相对位置误差');
grid on;ylim([-0.01,0.01]);

figure%理论相对位置与纯解算相对位置对比
subplot(311)
plot(T(2:4:end),R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,1),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('x/m');
title('天线中心相对位置R');
grid on;
ylim([0.48,0.52]);
subplot(312)
plot(T(2:4:end),R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,2),'--','LineWidth',2)
grid on;
legend('真实值','计算值');
xlabel('时间/s');ylabel('y/m');
title('天线中心相对位置R');ylim([-0.01,0.01]);
subplot(313)
plot(T(2:4:end),R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,3),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('z/m');
title('天线中心相对位置R');
grid on;