function avr=my_CellAvr(Cell)
% �������ƣ�my_CellAvr 
% �������ܣ�Ԫ���ֵ
% ���룺Cell      :Ԫ��
% �����avr       :��ֵ
[M,N]=size(Cell);
avr=zeros(size(Cell{1,1}));
for m=1:M
    for n=1:N
        avr=avr+Cell{m,n};
    end
end
avr=avr/(M*N);