function kfft=my_kfftA15(eth,Cmn,Csn,k,ts)
% �������ƣ�my_kfftA15
% �������ܣ���̬ƥ���˲�״̬ת�ƾ���
% ���룺eth      :���ߵ���ṹ�壬���Եõ�wnie,wnen,wnin
%       Csn      :�ӹߵ����ߵ���ϵ������任��3*3
%       fs       :�ӹ߱���3*1
%       ts       :�������
% �����kfft:״̬ת�ƾ����������ݾ���Ƚṹ��
%% �����һ��
A11=-askew(eth.wnin);A12=-Csn;A13_A15=zeros(3,9);
B11=-Csn;B12=zeros(3);
%% ����ڶ��С������С�������
A21_A3_5=zeros(6,15);A41_A4_4=zeros(3,12);A45=eye(3);
B21_B42=zeros(9,6);
%% ���������
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