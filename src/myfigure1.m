if flag.EnReli
    if flag.Ceil==0
        R=SINS1.R;%����λ��
        RR=SINSR1.Rall(1:2:end,:);
        T=ts:ts:(Sinf.len*ts);
    else
        R=SINS{flag.flagSINS_M,flag.flagSINS_N}.R;%����λ��
        RR=SINSR{flag.flagSINS_M,flag.flagSINS_N}.Rall(1:2:end,:);
        T=ts:ts:(SinfCell{1,1}.len*ts);
    end
    figure%�������λ���봿�������λ�öԱ�
    subplot(311)
    plot(T(2:4:end),R(1,2:4:end),'LineWidth',2)
    hold on;
    plot(T(2:4:end),RR(:,1),'--','LineWidth',2)
    legend('��ʵֵ','����ֵ');
    xlabel('ʱ��/s');ylabel('x/m');
    title('���ӹߵ����λ��R');
    grid on;
    subplot(312)
    plot(T(2:4:end),R(2,2:4:end),'LineWidth',2)
    hold on;
    plot(T(2:4:end),RR(:,2),'--','LineWidth',2)
    grid on;
    legend('��ʵֵ','����ֵ');
    xlabel('ʱ��/s');ylabel('y/m');
    title('���ӹߵ����λ��R');
    subplot(313)
    plot(T(2:4:end),R(3,2:4:end),'LineWidth',2)
    hold on;
    plot(T(2:4:end),RR(:,3),'--','LineWidth',2)
    legend('��ʵֵ','����ֵ');
    xlabel('ʱ��/s');ylabel('z/m');
    title('���ӹߵ����λ��R');
    grid on;
end