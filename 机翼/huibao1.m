close all
% load result_sim100
disp('max')
disp(max(result,[],2))
figure
subplot(131)
hist(result(1,:));
xlabel("Xλ������ֵ��/mm");ylabel("Ƶ��");
hold on;
subplot(132)
hist(result(2,:));
xlabel("Yλ������ֵ��/mm");ylabel("Ƶ��");
hold on;
subplot(133)
hist(result(3,:));
xlabel("Zλ������ֵ��/mm");ylabel("Ƶ��");
figure
subplot(131)
hist(result(4,:));
xlabel("X��̬����ֵ��/��");ylabel("Ƶ��");
hold on;
subplot(132)
hist(result(5,:));
xlabel("Y��̬����ֵ��/��");ylabel("Ƶ��");
hold on;
subplot(133)
hist(result(6,:));
xlabel("Z��̬����ֵ��/��");ylabel("Ƶ��");
%
figure
subplot(131)
hist(result(7,:));
xlabel("Xλ������ֵ��/mm");ylabel("Ƶ��");
hold on;
subplot(132)
hist(result(8,:));
xlabel("Yλ������ֵ��/mm");ylabel("Ƶ��");
hold on;
subplot(133)
hist(result(9,:));
xlabel("Zλ������ֵ��/mm");ylabel("Ƶ��");
figure
subplot(131)
hist(result(10,:));
xlabel("X��̬����ֵ��/��");ylabel("Ƶ��");
hold on;
subplot(132)
hist(result(11,:));
xlabel("Y��̬����ֵ��/��");ylabel("Ƶ��");
hold on;
subplot(133)
hist(result(12,:));
xlabel("Z��̬����ֵ��/��");ylabel("Ƶ��");
%
figure
subplot(131)
hist(result(13,:));
xlabel("Xλ�������/mm");ylabel("Ƶ��");
hold on;
subplot(132)
hist(result(14,:));
xlabel("Yλ�������/mm");ylabel("Ƶ��");
hold on;
subplot(133)
hist(result(15,:));
xlabel("Zλ�������/mm");ylabel("Ƶ��");
figure
subplot(131)
hist(result(16,:));
xlabel("X��̬�����/��");ylabel("Ƶ��");
hold on;
subplot(132)
hist(result(17,:));
xlabel("Y��̬�����/��");ylabel("Ƶ��");
hold on;
subplot(133)
hist(result(18,:));
xlabel("Z��̬�����/��");ylabel("Ƶ��");
Vmean=mean(abs(result),2)
VErr=my_getVErr(abs(result),1,size(result,2),2)