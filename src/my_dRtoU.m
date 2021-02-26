function U=my_dRtoU(R,varargin)
% 函数名称：my_dRtoU
% 函数功能：反解算相对位置矢量微分方程，对应my_UtodiffU
% 输入：R      :当前导航系下相对位置
%       varargin:函数除微分项其余参数,以元组给入          
%              每一个参数中依次存放dR 位置微分
%                                  wim 主惯相对于惯性坐标系角速度在主惯坐标系的投影
% 输出：dR:当前导航系下相对位置微分
dR=varargin{1};
wim=varargin{2};
U=dR+cross(wim,R);