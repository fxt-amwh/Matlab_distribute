function kfft=my_kfftA15(eth,Cmn,Csn,k,ts)
% 函数名称：my_kfftA15
% 函数功能：姿态匹配滤波状态转移矩阵
% 输入：eth      :主惯地球结构体，可以得到wnie,wnen,wnin
%       Csn      :子惯到主惯导航系的坐标变换阵3*3
%       fs       :子惯比力3*1
%       ts       :采样间隔
% 输出：kfft:状态转移矩阵、噪声传递矩阵等结构体
%% 矩阵第一行
A11=-askew(eth.wnin);A12=-Csn;A13_A15=zeros(3,9);
B11=-Csn;B12=zeros(3);
%% 矩阵第二行、第三行、第四行
A21_A3_5=zeros(6,15);A41_A4_4=zeros(3,12);A45=eye(3);
B21_B42=zeros(9,6);
%% 矩阵第五行
A51_A53=zeros(3,9);A54=-diag(k.^2);A55=-2*diag(k);
B51=zeros(3);B52=eye(3);
%% 
kfft.A=[ A11,A12,A13_A15;
    A21_A3_5;
    A41_A4_4,A45;
    A51_A53,A54,A55];
kfft.phi=eye(15)+kfft.A*ts;
kfft.B=[B11,B12;
    B21_B42;
    B51,B52];
kfft.Gammak=(eye(15)+kfft.A*ts/2)*kfft.B;
kfft.H=[-eye(3),zeros(3,3),Cmn,Cmn,zeros(3,3)];