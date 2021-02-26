function [qnb, vn, pos, eth] = my_SINSVgetResult(qnb, vn, pos, wm, vm, ts)
% 函数名称：my_SINSVgetResult
% 函数功能：由子惯输出进行相对导航
% 输入：MINS     ：主惯结构体
%       SINS     ：子惯结构体
%       varargin:函数除微分项其余参数,以元组给入          
%              每一个参数中依次存放Cms，（主到子坐标变换阵，与Csm区分）3*3
%                                  dU  速度微分
%                                  fm  主惯比力
%                                  wim 主惯相对于惯性坐标系速度在主惯坐标系的投影
% 输出：INSSR:相对导航结果
nn = size(wm,1);  nts = nn*ts;
[phim, dvbm] = cnscl(wm, vm);  % 圆锥误差/划船误差补偿
eth = earth(pos, vn);  % 地球相关参数计算
vn1 = vn + rv2m(-eth.wnin*nts/2)*qmulv(qnb,dvbm);  % 速度更新
vn = (vn+vn1)/2;
pos = pos + [vn(2)/eth.RMh;vn(1)/eth.clRNh;vn(3)]*nts;  vn = vn1;  % 位置更新
qnb = qmul(rv2q(-eth.wnin*nts), qmul(qnb, rv2q(phim)));  % 姿态更新
qnb = qnormlz(qnb);