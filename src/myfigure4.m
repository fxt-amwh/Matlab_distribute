if flag.Ceil==0
    atttrue=SINS1.atttrue(:,2:4:end)/arcdeg;%����λ��
    atttrueFR=SINSFR1.attall(:,1:2:end)/arcdeg;
    T=ts:ts:(Sinf.len*ts);
else
    atttrue=SINS{flag.flagSINS_M,flag.flagSINS_N}.atttrue(:,2:4:end)/arcdeg;%����λ��
    atttrueFR=SINSFR{flag.flagSINS_M,flag.flagSINS_N}.attall(:,1:2:end)/arcdeg;
    T=ts:ts:(SinfCell{1,1}.len*ts);
end
figure%���۽Ƕ����˲�����ǶȶԱ�
subplot(311)
plot(T(2:4:end-1),atttrue(1,:),'LineWidth',2)
hold on;
plot(T(2:4:end-1),atttrueFR(1,:),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('������/��');
title('���ӹߵ������̬���˲���');

grid on;
subplot(312)
plot(T(2:4:end-1),atttrue(2,:),'LineWidth',2)
hold on;
plot(T(2:4:end-1),atttrueFR(2,:),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('��ת��/��');
title('���ӹߵ������̬���˲���');

grid on;
subplot(313)
plot(T(2:4:end-1),atttrue(3,:),'LineWidth',2)
hold on;
plot(T(2:4:end-1),atttrueFR(3,:),'--','LineWidth',2)
legend('��ʵֵ','����ֵ');
xlabel('ʱ��/s');ylabel('�����/��');
title('���ӹߵ������̬���˲���');

grid on;