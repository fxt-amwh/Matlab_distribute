function kfft=my_kfftV12(eth,Csn,fs,ts)
% �������ƣ�my_kfftV12
% �������ܣ��ٶ�ƥ���˲�״̬ת�ƾ���
% ���룺eth      :���ߵ���ṹ�壬���Եõ�wnie,wnen,wnin
%       Csn      :�ӹߵ����ߵ���ϵ������任��3*3
%       fs       :�ӹ߱���3*1
%       ts       :�������
% �����kfft:״̬ת�ƾ����������ݾ���Ƚṹ��
%% �����һ��
A11=-askew(2*eth.wnie+eth.wnen);A12=askew(Csn*fs);A13=Csn;A14=zeros(3);
B11=Csn;B12=zeros(3);
%% ����ڶ���
A21=zeros(3);A22=askew(-eth.wnin);A23=zeros(3);A24=-Csn;
B21=zeros(3);B22=-Csn;
%% ��������С�������
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