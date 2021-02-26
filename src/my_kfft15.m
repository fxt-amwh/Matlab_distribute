function kfft=my_kfft15(wim,Csm,fs,ts)
% �������ƣ�my_kfft15
% �������ܣ���Ե����˲�״̬ת�ƾ���
% ���룺wim      :���߽��ٶ�3*1
%       Csm      :�ӹߵ����ߵ�����任��3*3
%       fs       :�ӹ߱���3*1
%       ts       :�������
% �����kfft:״̬ת�ƾ����������ݾ���Ƚṹ��
%% �����һ��
A11=-askew(wim);A12=zeros(3);A13=zeros(3);A14=-Csm;A15=zeros(3);
B11=-Csm;B12=zeros(3);
%% ����ڶ���
A21=askew(Csm*fs);A22=-askew(wim);A23=zeros(3);A24=zeros(3);A25=Csm;
B21=zeros(3);B22=Csm;
%% ���������
A31=zeros(3);A32=eye(3);A33=-askew(wim);A34=zeros(3);A35=zeros(3);
B31_B52=zeros(9,6);
%% ��������С�������
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