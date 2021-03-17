function err=my_getSErrBE(ex,beg,End,dir)
% 函数名称：my_getSErr
% 函数功能：计算最终的误差
% 输入：ex      :差值
%       dir     :方向，1：按列计算，2按行计算
% 输出：INSSFR:相对导航结果
len=End-beg+1;
if dir==1
    err=sum(abs(ex(beg:End,:)),dir)/len;
elseif dir==2
    err=sum(abs(ex(:,beg:End)),dir)/len;
end