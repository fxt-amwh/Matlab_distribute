%% 对比1m处子惯与2m处子惯原始数据
figure
subplot(2,1,1)
plot(SINS{1, 1}.U','DisplayName','SINS{1, 1}.U')
ylim([-0.02,0.04])
subplot(2,1,2)
plot(SINS{1, 4}.U','DisplayName','SINS{1, 1}.U')
ylim([-0.02,0.04])
figure
subplot(2,1,1)
plot(2*SINS{1, 1}.U','DisplayName','SINS{1, 1}.U')
ylim([-0.02,0.04])
subplot(2,1,2)
plot(SINS{1, 4}.U','DisplayName','SINS{1, 1}.U')
ylim([-0.02,0.04])
figure
subplot(3,1,1)
plot(2*SINS{1, 1}.U(1,:),'DisplayName','SINS{1, 1}.U')
hold on;
plot(SINS{1, 4}.U(1,:),'DisplayName','SINS{1, 1}.U')
subplot(3,1,2)
plot(2*SINS{1, 1}.U(2,:),'DisplayName','SINS{1, 1}.U')
hold on;
plot(SINS{1, 4}.U(2,:),'DisplayName','SINS{1, 1}.U')
subplot(3,1,3)
plot(2*SINS{1, 1}.U(3,:),'DisplayName','SINS{1, 1}.U')
hold on;
plot(SINS{1, 4}.U(3,:),'DisplayName','SINS{1, 1}.U')
%% 对比1m处子惯与2m处子惯原始数据 全部转换到2m处
 klist=R_lw(1,1)./Sinf.Rlist(1,:);
k=2;%wms系数
eb=SERR{1,1}.eb;
web=SERR{1,1}.web;
sts = sqrt(ts);
wm = wm + [ eb(1) + sts*web(1)*randn(1,1)/ts,eb(2) + sts*web(2)*randn(1,1)/ts,eb(3) + sts*web(3)*randn(1,1)/ts ];
% fm = fm + [ db(1) + sts*wdb(1)*randn(m,1)/ts,db(2) + sts*wdb(2)*randn(m,1)/ts,db(3) + sts*wdb(3)*randn(m,1)/ts ];
wms_m1=zeros(length(SINS{1, 1}.wis(1,:)),3);
wms_m2_true=zeros(length(SINS{1, 1}.wis(1,:)),3);
wms_m1_3=zeros(length(SINS{1, 1}.wis(1,:)),3);
wms_m1_7=zeros(length(SINS{1, 1}.wis(1,:)),3);
wms_m2=zeros(length(SINS{1, 1}.wis(1,:)),3);
wms_mave=zeros(length(SINS{1, 1}.wis(1,:)),3);
wms_m2_noerr=zeros(length(SINS{1, 1}.wis(1,:)),3);
lendata=length(SINS{1, 1}.wis(1,:))-1;
for i=1:length(SINS{1, 1}.wis(1,:))
    fitatterr=0*[0 1 1 1
                1 0 1 0
                1 1 0 0 ]*d2g;
    errC1=a2mat(fitatterr(:,1)); errC2=a2mat(fitatterr(:,2)); errC3=a2mat(fitatterr(:,3)); errC4=a2mat(fitatterr(:,4)); 
    wis2_true=SINS{1, 4}.wis(:,i);
    wms_m2_true(i,:)=SINS{1, 4}.Cms{i}'*wis2_true-MINS.wim(:,i);
    wis1=SINS{1, 1}.wis(:,i)+[ eb(1) + sts*web(1)*randn(1,1)/ts;eb(2) + sts*web(2)*randn(1,1)/ts;eb(3) + sts*web(3)*randn(1,1)/ts ];
    wis1_3=SINS{1, 2}.wis(:,i)+[ eb(1) + sts*web(1)*randn(1,1)/ts;eb(2) + sts*web(2)*randn(1,1)/ts;eb(3) + sts*web(3)*randn(1,1)/ts ];
    wis1_7=SINS{1, 3}.wis(:,i)+[ eb(1) + sts*web(1)*randn(1,1)/ts;eb(2) + sts*web(2)*randn(1,1)/ts;eb(3) + sts*web(3)*randn(1,1)/ts ];
    wis2=SINS{1, 4}.wis(:,i)+[ eb(1) + sts*web(1)*randn(1,1)/ts;eb(2) + sts*web(2)*randn(1,1)/ts;eb(3) + sts*web(3)*randn(1,1)/ts ];
    wms_m1(i,:)=SINS{1, 1}.Cms{i}'*errC1*wis1-MINS.wim(:,i);
    wms_m1_3(i,:)=SINS{1, 2}.Cms{i}'*errC2*wis1_3-MINS.wim(:,i);
    wms_m1_7(i,:)=SINS{1, 3}.Cms{i}'*errC3*wis1_7-MINS.wim(:,i);
    wms_m2(i,:)=SINS{1, 4}.Cms{i}'*errC4*wis2-MINS.wim(:,i);
    wms_m2_noerr(i,:)=SINS{1, 4}.Cms{i}'*wis2-MINS.wim(:,i);
    wms_mave(i,:)=(wms_m1(i,:).*[1,klist(1,1),1]+wms_m1_3(i,:).*[1,klist(1,2),1]+wms_m1_7(i,:).*[1,klist(1,3),1]+wms_m2(i,:).*[1,klist(1,4),1])/4;
end
figure
subplot(3,1,1)
plot(wms_m1(1:lendata,1)/arcdeg,'DisplayName','SINS{1, 1}.U')
hold on;
plot(wms_m2(1:lendata,1)/arcdeg,'DisplayName','SINS{1, 1}.U')
hold on;
plot(wms_mave(1:lendata,1)/arcdeg,'DisplayName','SINS{1, 1}.U')
subplot(3,1,2)
plot(k*wms_m1(1:lendata,2)/arcdeg,'DisplayName','SINS{1, 1}.U','LineWidth',2)
hold on;
plot(wms_m2(1:lendata,2)/arcdeg,'--','DisplayName','SINS{1, 1}.U','LineWidth',2)
hold on;
plot(wms_mave(1:lendata,2)/arcdeg,'.-','DisplayName','SINS{1, 1}.U','LineWidth',1)
subplot(3,1,3)
plot(wms_m1(1:lendata,3)/arcdeg,'DisplayName','SINS{1, 1}.U')
hold on;
plot(wms_m2(1:lendata,3)/arcdeg,'DisplayName','SINS{1, 1}.U')
hold on;
plot(wms_mave(1:lendata,3)/arcdeg,'DisplayName','SINS{1, 1}.U')
figure
subplot(3,1,1)
plot((wms_m2_noerr(1:lendata,1)-wms_m2_true(1:lendata,1))/arcdeg,'DisplayName','SINS{1, 1}.U')
hold on;
plot((wms_mave(1:lendata,1)-wms_m2_true(1:lendata,1))/arcdeg,'DisplayName','SINS{1, 1}.U')
legend("无初始误差时噪声数据","有安装误差时多节点融合数据误差")
subplot(3,1,2)
plot((wms_m2_noerr(1:lendata,2)-wms_m2_true(1:lendata,2))/arcdeg,'DisplayName','SINS{1, 1}.U')
hold on;
plot((wms_mave(1:lendata,2)-wms_m2_true(1:lendata,2))/arcdeg,'DisplayName','SINS{1, 1}.U')
legend("无初始误差时噪声数据","有安装误差时多节点融合数据误差")
subplot(3,1,3)
plot((wms_m2_noerr(1:lendata,3)-wms_m2_true(1:lendata,3))/arcdeg,'DisplayName','SINS{1, 1}.U')
hold on;
plot((wms_mave(1:lendata,3)-wms_m2_true(1:lendata,3))/arcdeg,'DisplayName','SINS{1, 1}.U')
legend("无初始误差时噪声数据","有安装误差时多节点融合数据误差")
var(wms_m1(1:lendata,3)/arcdeg)
var(wms_mave(1:lendata,3)/arcdeg)
%% 对比1m处子惯与2m处子惯原始数据
figure
subplot(2,1,1)
plot(SINS{1, 1}.fs','DisplayName','SINS{1, 1}.U')
subplot(2,1,2)
plot(SINS{1, 4}.fs','DisplayName','SINS{1, 1}.U')