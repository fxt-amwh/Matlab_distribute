close all;
figure
plot(NSINSFR.Rloop_k(2,:)*180/pi);hold on;
plot(NSINSFR.Rloop_k(5,:)*180/pi);hold on;
plot(NSINSFR.Rloop_k(8,:)*180/pi);hold on;
plot(NSINSFR.Rloop_k(11,:)*180/pi);
legend("1","2","3","4");
title("分别挠曲角");
figure
plot(NSINSFR.uf(2,:)*180/pi)
title("融合挠曲角");
figure
plot(NSINSFR.Rloopk(1,:));hold on;
plot(NSINSFR.Rloopk(4,:));hold on;
plot(NSINSFR.Rloopk(7,:));hold on;
plot(NSINSFR.Rloopk(10,:));
legend("1","2","3","4");
title("分别挠曲位移x");
figure
plot(NSINSFR.Rloopk(2,:));hold on;
plot(NSINSFR.Rloopk(5,:));hold on;
plot(NSINSFR.Rloopk(8,:));hold on;
plot(NSINSFR.Rloopk(11,:));
legend("1","2","3","4");
title("分别挠曲位移y");
figure
plot(NSINSFR.Rloopk(3,:));hold on;
plot(NSINSFR.Rloopk(6,:));hold on;
plot(NSINSFR.Rloopk(9,:));hold on;
plot(NSINSFR.Rloopk(12,:));
legend("1","2","3","4");
title("分别挠曲位移z");