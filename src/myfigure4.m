if flag.Ceil==0
    atttrue=SINS1.atttrue(:,2:4:end)/arcdeg;%理论位置
    atttrueFR=SINSFR1.attall(:,1:2:end)/arcdeg;
    T=ts:ts:(Sinf.len*ts);
else
    atttrue=SINS{flag.flagSINS_M,flag.flagSINS_N}.atttrue(:,2:4:end)/arcdeg;%理论位置
    atttrueFR=SINSFR{flag.flagSINS_M,flag.flagSINS_N}.attall(:,1:2:end)/arcdeg;
    T=ts:ts:(SinfCell{1,1}.len*ts);
end
figure%理论角度与滤波解算角度对比
subplot(311)
plot(T(2:4:end-1),atttrue(1,:),'LineWidth',2)
hold on;
plot(T(2:4:end-1),atttrueFR(1,:),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('俯仰角/°');
title('主子惯导相对姿态（滤波）');

grid on;
subplot(312)
plot(T(2:4:end-1),atttrue(2,:),'LineWidth',2)
hold on;
plot(T(2:4:end-1),atttrueFR(2,:),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('滚转角/°');
title('主子惯导相对姿态（滤波）');

grid on;
subplot(313)
plot(T(2:4:end-1),atttrue(3,:),'LineWidth',2)
hold on;
plot(T(2:4:end-1),atttrueFR(3,:),'--','LineWidth',2)
legend('真实值','计算值');
xlabel('时间/s');ylabel('航向角/°');
title('主子惯导相对姿态（滤波）');

grid on;