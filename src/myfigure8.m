if flag.Ceil==0
    atttrue=SINS1.atttrue(:,2:4:end)/arcdeg;%����λ��
    atttrueFR=NSINSFR.attall(:,1:2:end)/arcdeg;
    T=ts:ts:(Sinf.len*ts);
else
    atttrue=SINS{flag.flagSINS_M,flag.flagSINS_N}.atttrue(:,2:4:end)/arcdeg;%����λ��
    atttrueFR=NSINSFR.attall(:,1:2:end)/arcdeg;
    T=ts:ts:(SinfCell{1,1}.len*ts);
end
figure
subplot(311)
plot(T(2:4:end-1),atttrueFR(1,:)-atttrue(1,:),'LineWidth',2)
xlabel('ʱ��/s');ylabel('������/��');title('���ӹߵ������̬���ںϣ�');grid on;ylim([-0.1,0.1]);
subplot(312)
plot(T(2:4:end-1),atttrueFR(2,:)-atttrue(2,:),'LineWidth',2)
xlabel('ʱ��/s');ylabel('��ת��/��');title('���ӹߵ������̬���ںϣ�');grid on;ylim([-0.1,0.1]);
subplot(313)
plot(T(2:4:end-1),atttrueFR(3,:)-atttrue(3,:),'LineWidth',2)
xlabel('ʱ��/s');ylabel('�����/��');title('���ӹߵ������̬���ںϣ�');grid on;ylim([-0.1,0.1]);

figure%���۽Ƕ��봿��Խ���ǶȶԱ�
subplot(311)
plot(T(2:4:end),atttrue(1,:),'LineWidth',2)
hold on;
plot(T(2:4:end),atttrueFR(1,:),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('������/��');
title('���ӹߵ������̬');
grid on;ylim([0,2]);

subplot(312)
plot(T(2:4:end),atttrue(2,:),'LineWidth',2)
hold on;
plot(T(2:4:end),atttrueFR(2,:),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('��ת��/��');
title('���ӹߵ������̬');
grid on;

subplot(313)
plot(T(2:4:end),atttrue(3,:),'LineWidth',2)
hold on;
plot(T(2:4:end),atttrueFR(3,:),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('�����/��');
title('���ӹߵ������̬');
grid on;ylim([2,4]);