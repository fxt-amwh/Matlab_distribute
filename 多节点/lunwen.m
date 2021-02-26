if flag.Ceil==0
    R=SINS1.R;%����λ��
    RFR=NSINSFR.Rall(1:2:end,:);
    T=ts:ts:(Sinf.len*ts);
else
    R=SINS{flag.flagSINS_M,flag.flagSINS_N}.R;%����λ��
    RFR=NSINSFR.Rall(1:2:end,:);
    T=ts:ts:(SinfCell{1,1}.len*ts);
end
figure
subplot(311)
plot(T(2:4:end),RFR(:,1)'-R(1,2:4:end),'LineWidth',2)
legend('X');
xlabel('ʱ��/s');ylabel('x/m');
title('�����������λ�����');
grid on;ylim([-0.01,0.01]);
subplot(312)
plot(T(2:4:end),RFR(:,2)'-R(2,2:4:end),'LineWidth',2)
grid on;
legend('Y');
xlabel('ʱ��/s');ylabel('y/m');
title('�����������λ�����');ylim([-0.01,0.01]);
subplot(313)
plot(T(2:4:end),RFR(:,3)'-R(3,2:4:end),'LineWidth',2)
legend('Z');
xlabel('ʱ��/s');ylabel('z/m');
title('�����������λ�����');
grid on;ylim([-0.01,0.01]);

figure%�������λ���봿�������λ�öԱ�
subplot(311)
plot(T(2:4:end),R(1,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,1),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('x/m');
title('�����������λ��R');
grid on;
ylim([0.48,0.52]);
subplot(312)
plot(T(2:4:end),R(2,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,2),'--','LineWidth',2)
grid on;
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('y/m');
title('�����������λ��R');ylim([-0.01,0.01]);
subplot(313)
plot(T(2:4:end),R(3,2:4:end),'LineWidth',2)
hold on;
plot(T(2:4:end),RFR(:,3),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('z/m');
title('�����������λ��R');
grid on;