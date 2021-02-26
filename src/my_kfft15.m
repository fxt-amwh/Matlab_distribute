function kfft=my_kfft15(wim,Csm,fs,ts)
% 函数名称：my_kfft15
% 函数功能：相对导航滤波状态转移矩阵
% 输入：wim      :主惯角速度3*1
%       Csm      :子惯到主惯的坐标变换阵3*3
%       fs       :子惯比力3*1
%       ts       :采样间隔
% 输出：kfft:状态转移矩阵、噪声传递矩阵等结构体
%% 矩阵第一行
A11=-askew(wim);A12=zeros(3);A13=zeros(3);A14=-Csm;A15=zeros(3);
B11=-Csm;B12=zeros(3);
%% 矩阵第二行
A21=askew(Csm*fs);A22=-askew(wim);A23=zeros(3);A24=zeros(3);A25=Csm;
B21=zeros(3);B22=Csm;
%% 矩阵第三行
A31=zeros(3);A32=eye(3);A33=-askew(wim);A34=zeros(3);A35=zeros(3);
B31_B52=zeros(9,6);
%% 矩阵第四行、第五行
A41_A55=zeros(6,15);
%% 
kfft.A=[ A11,A12,A13,A14,A15;
    A21,A22,A23,A24,A25;
    A31,A32,A33,A34,A35
    A41_A55];
kfft.phi=eye(15)+kfft.A*ts;
kfft.B=[B11,B12;
    B21,B22;
    B31_B52];
kfft.Gammak=(eye(15)+kfft.A*ts/2)*kfft.B;
kfft.H=[zeros(3,6),eye(3),zeros(3,6)];