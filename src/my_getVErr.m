function err=my_getVErr(ex,beg,End,dir)
% 函数名称：my_getVErr
% 函数功能：计算最终的均方差
% 输入：beg      :差值
%       end     :方向，1：按列计算，2按行计算
%       dir     :方向，1：按列计算，2按行计算
% 输出：INSSFR:相对导航结果
if dir==1
    err=sqrt(var(ex(beg:End,:),0,dir));
elseif dir==2
    err=sqrt(var(ex(:,beg:End),0,dir));
end