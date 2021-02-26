function fs=my_dUtofs(U,varargin)
% 函数名称：my_dUtofs
% 函数功能：反解算相对导航速度微分方程,与my_UtodiffU对应
% 输入：U      :当前导航系下相对速度
%       varargin:函数除微分项其余参数,以元组给入          
%              每一个参数中依次存放Cms，（主到子坐标变换阵，与Csm区分）3*3
%                                  dU  速度微分
%                                  fm  主惯比力
%                                  wim 主惯相对于惯性坐标系速度在主惯坐标系的投影
% 输出：dU:当前导航系下相对速度微分
Cms=varargin{1};
dU=varargin{2};
fm=varargin{3};
wim=varargin{4};
fs=Cms*(dU+fm+cross(wim,U));