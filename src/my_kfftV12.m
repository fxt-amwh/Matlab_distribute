function kfft=my_kfftV12(eth,Csn,fs,ts)
% 函数名称：my_kfftV12
% 函数功能：速度匹配滤波状态转移矩阵
% 输入：eth      :主惯地球结构体，可以得到wnie,wnen,wnin
%       Csn      :子惯到主惯导航系的坐标变换阵3*3
%       fs       :子惯比力3*1
%       ts       :采样间隔
% 输出：kfft:状态转移矩阵、噪声传递矩阵等结构体
%% 矩阵第一行
A11=-askew(2*eth.wnie+eth.wnen);A12=askew(Csn*fs);A13=Csn;A14=zeros(3);
B11=Csn;B12=zeros(3);
%% 矩阵第二行
A21=zeros(3);A22=askew(-eth.wnin);A23=zeros(3);A24=-Csn;
B21=zeros(3);B22=-Csn;
%% 矩阵第三行、第四行
A31_A44=zeros(6,12);
B31_B42=zeros(6,6);
%% 
kfft.A=[ A11,A12,A13,A14;
    A21,A22,A23,A24;
    A31_A44];
kfft.phi=eye(12)+kfft.A*ts;
kfft.B=[B11,B12;
    B21,B22;
    B31_B42];
kfft.Gammak=(eye(12)+kfft.A*ts/2)*kfft.B;
kfft.H=[eye(3),zeros(3,9)];