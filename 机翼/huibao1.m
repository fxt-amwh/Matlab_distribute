close all
% load result_sim100
disp('max')
disp(max(result,[],2))
figure
subplot(131)
hist(result(1,:));
xlabel("X位置误差（终值）/mm");ylabel("频数");
hold on;
subplot(132)
hist(result(2,:));
xlabel("Y位置误差（终值）/mm");ylabel("频数");
hold on;
subplot(133)
hist(result(3,:));
xlabel("Z位置误差（终值）/mm");ylabel("频数");
figure
subplot(131)
hist(result(4,:));
xlabel("X姿态误差（终值）/分");ylabel("频数");
hold on;
subplot(132)
hist(result(5,:));
xlabel("Y姿态误差（终值）/分");ylabel("频数");
hold on;
subplot(133)
hist(result(6,:));
xlabel("Z姿态误差（终值）/分");ylabel("频数");
%
figure
subplot(131)
hist(result(7,:));
xlabel("X位置误差（均值）/mm");ylabel("频数");
hold on;
subplot(132)
hist(result(8,:));
xlabel("Y位置误差（均值）/mm");ylabel("频数");
hold on;
subplot(133)
hist(result(9,:));
xlabel("Z位置误差（均值）/mm");ylabel("频数");
figure
subplot(131)
hist(result(10,:));
xlabel("X姿态误差（均值）/分");ylabel("频数");
hold on;
subplot(132)
hist(result(11,:));
xlabel("Y姿态误差（均值）/分");ylabel("频数");
hold on;
subplot(133)
hist(result(12,:));
xlabel("Z姿态误差（均值）/分");ylabel("频数");
%
figure
subplot(131)
hist(result(13,:));
xlabel("X位置误差（方差）/mm");ylabel("频数");
hold on;
subplot(132)
hist(result(14,:));
xlabel("Y位置误差（方差）/mm");ylabel("频数");
hold on;
subplot(133)
hist(result(15,:));
xlabel("Z位置误差（方差）/mm");ylabel("频数");
figure
subplot(131)
hist(result(16,:));
xlabel("X姿态误差（方差）/分");ylabel("频数");
hold on;
subplot(132)
hist(result(17,:));
xlabel("Y姿态误差（方差）/分");ylabel("频数");
hold on;
subplot(133)
hist(result(18,:));
xlabel("Z姿态误差（方差）/分");ylabel("频数");
Vmean=mean(abs(result),2)
VErr=my_getVErr(abs(result),1,size(result,2),2)