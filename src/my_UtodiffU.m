function dU=my_UtodiffU(U,u)
% 函数名称：my_UtodiffU
% 函数功能：相对导航速度微分方程
% 输入：U      :当前导航系下相对速度
%       u      :函数除微分项其余参数,以元组给入          
%              每一个参数中依次存放Cms，（主到子坐标变换阵，与Csm区分）3*3
%                                  fs  子惯比力
%                                  fm  主惯比力
%                                  wim 主惯相对于惯性坐标系速度在主惯坐标系的投影
% 输出：dU:当前导航系下相对速度微分
Cms=u{1};
fs=u{2};
fm=u{3};
wim=u{4};
Csm=Cms';
dU=Csm*fs-fm-cross(wim,U);
% dU=Csm*fs-fm;
