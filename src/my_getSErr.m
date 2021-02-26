function err=my_getSErr(ex,len,dir)
% 函数名称：my_getSErr
% 函数功能：计算最终的误差
% 输入：ex      :差值
%       dir     :方向，1：按列计算，2按行计算
% 输出：INSSFR:相对导航结果
len=floor(len);
if dir==1
    err=sum(ex(end-len+1:end,:),dir)/len;
elseif dir==2
    err=sum(ex(:,end-len+1:end),dir)/len;
end