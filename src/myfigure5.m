if flag.Ceil==0
    R=SINS1.R;%����λ��
    RFR=SINSFR1.Rall(1:2:end,:);
    T=ts:ts:(Sinf.len*ts);
else
    R=SINS{flag.flagSINS_M,flag.flagSINS_N}.R;%����λ��
    RFR=SINSFR{flag.flagSINS_M,flag.flagSINS_N}.Rall(1:2:end,:);
    T=ts:ts:(SinfCell{1,1}.len*ts);
end
figure
subplot(311)
plot(T(2:4:end),RFR(:,1)'-R(1,2:4:end),'LineWidth',2)
legend('X');
xlabel('ʱ��/s');ylabel('x/m');
title('���ӹߵ����λ�����');
ylim([-0.001,0.001]);
grid on;
subplot(312)
plot(T(2:4:end),RFR(:,2)'-R(2,2:4:end),'LineWidth',2)
grid on;
legend('Y');
xlabel('ʱ��/s');ylabel('y/m');
title('���ӹߵ����λ�����');
ylim([-0.001,0.001]);
subplot(313)
plot(T(2:4:end),RFR(:,3)'-R(3,2:4:end),'LineWidth',2)
legend('Z');
xlabel('ʱ��/s');ylabel('z/m');
title('���ӹߵ����λ�����');
ylim([-0.001,0.001]);
grid on;