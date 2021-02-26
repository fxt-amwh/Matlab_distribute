function [qnm,vnm,posm,qns,vns,poss,VL,ethm] = my_Vrelinsupdate(qnm,vnm,posm,wm,vm,qns,vns,poss,ws,vs,R,ts)
% 函数名称：my_Vrelinsupdate
% 函数功能：由子惯输出进行相对导航
% 输入：qnm     ：主惯姿态四元数 注意含义，由于历史原因，这里的qnm表示由主惯坐标系到导航坐标系四元数
%       vnm     ：主惯速度输出 
%       posm    ：主惯经纬高
%       wm      ：主惯陀螺增量
%       vm      ：主惯速度增量
%       qns     ：子惯姿态四元数 注意含义，由于历史原因，这里的qns表示由子惯坐标系到导航坐标系四元数
%       vns     ：子惯速度输出
%       poss    ：子惯经纬高
%       ws      ：子惯陀螺增量
%       vs      ：子惯速度增量
%       R       ：主子间杆臂
%       ts      ：采样时间
% 输出：qnm     ：主惯姿态四元数 注意含义，由于历史原因，这里的qnm表示由主惯坐标系到导航坐标系四元数
%       vnm     ：主惯速度输出
%       posm    ：主惯经纬高
%       qns     ：子惯姿态四元数 注意含义，由于历史原因，这里的qns表示由子惯坐标系到导航坐标系四元数
%       vns     ：子惯速度输出
%       poss    ：子惯经纬高
%       VL      ：主惯坐标系下子惯杆臂速度
%       ethm    ：主惯标准下地球相关参数
%主惯
nn = size(wm,1);  nts = nn*ts;
[phim, dvbm] = cnscl(wm, vm);  % 圆锥误差/划船误差补偿
ethm = earth(posm, vnm);  % 地球相关参数计算
vn1 = vnm + rv2m(-ethm.wnin*nts/2)*qmulv(qnm,dvbm);  % 速度更新
vnm = (vnm+vn1)/2;
posm = posm + [vnm(2)/ethm.RMh;vnm(1)/ethm.clRNh;vnm(3)]*nts;  vnm = vn1;  % 位置更新
qnm = qmul(rv2q(-ethm.wnin*nts), qmul(qnm, rv2q(phim)));  % 姿态更新
qnm = qnormlz(qnm);
%子惯
[phis, dvbs] = cnscl(ws, vs);  % 圆锥误差/划船误差补偿
eths = earth(poss, vns);  % 地球相关参数计算
vns1 = vns + rv2m(-eths.wnin*nts/2)*qmulv(qns,dvbs);  % 速度更新
vns = (vns+vns1)/2;
poss = poss + [vns(2)/eths.RMh;vns(1)/eths.clRNh;vns(3)]*nts;  vns = vns1;  % 位置更新
qns = qmul(rv2q(-eths.wnin*nts), qmul(qns, rv2q(phis)));  % 姿态更新
qns = qnormlz(qns);
%杆臂速度
w_MiM=mean(wm,1)/ts;%角增量转换为角速度
w_nie=ethm.wnie;
CnM=q2mat(qnm);
% VL=CMn*cross(w_MiM,R)'-cross(w_nie,CMn*R);
% VL=CMn*(askew(w_MiM'-CnM*w_nie))*R;
VL=CnM*cross(w_MiM,R)'-askew(w_nie)*CnM*R;