function KFinitCell=my_KFinitN(SERR,SINS,nts)
% �������ƣ�KFinitN
% �������ܣ����󿨶����˲���
% ���룺SERR      :�ӹ����Ԫ��
%       SINS      :�ӹ���ϢԪ��  
%       nts = 2*ts;�������Ͳ���ʱ��
% �����KFinitCell:�˲���Ԫ��
[M,N]=size(SERR);%Ԥ����������SINS�ӿ� M=1Ϊ��SINS��
KFinitCell=cell(M,N);
arcdeg = pi/180;
for m=1:M
    for n=1:N
        KFinit.Qk = diag([SERR{m,n}.web; SERR{m,n}.wdb;])^2*nts;
        KFinit.rk = [0.001;0.001;0.001]*sqrt(SINS{m,n}.R0'*SINS{m,n}.R0);  
        KFinit.Rk = diag(KFinit.rk)^2;
        KFinit.P0 = diag([[0.1;0.1;0.1]*arcdeg; [20;20;20]; [0.002;0.002;0.002]*sqrt(SINS{m,n}.R0'*SINS{m,n}.R0);
        SERR{m,n}.eb; SERR{m,n}.db])^2;
        KFinitCell{m,n}=KFinit;
    end
end