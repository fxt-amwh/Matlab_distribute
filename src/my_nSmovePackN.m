function SmoveCell=my_nSmovePackN(Sinf)
% 函数名称：my_nSmovePackN
% 函数功能：反解算相对导航子惯无噪声数据解封 面阵
% 输入：Sinf      :子惯信息配置矩阵 元组
% 输出：SmoveCell :子惯结构体元组
M=size(Sinf,2);%预留升级面阵SINS接口 M=1为线SINS组
N=size(Sinf{1}.Rlist,2);
SmoveCell=cell(M,N);
for m=1:M
    SmoveCellAllN=my_nSmovePack(Sinf{m});
    for n=1:N
        SmoveCell{m,n}=SmoveCellAllN{1,n};
    end
end