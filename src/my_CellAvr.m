function avr=my_CellAvr(Cell)
% 函数名称：my_CellAvr 
% 函数功能：元组均值
% 输入：Cell      :元组
% 输出：avr       :均值
[M,N]=size(Cell);
avr=zeros(size(Cell{1,1}));
for m=1:M
    for n=1:N
        avr=avr+Cell{m,n};
    end
end
avr=avr/(M*N);