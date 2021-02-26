function SmoveCell=my_nSmovePack(Sinf)
% 函数名称：my_nSmovePack
% 函数功能：反解算相对导航子惯无噪声数据解封
% 输入：Sinf      :子惯信息配置矩阵
% 输出：SmoveCell :子惯结构体元组
arcdeg = pi/180;
n=size(Sinf.Rlist,2);
SmoveCell=cell(1,n);
for loopi=1:n
    Smove.um=Sinf.ulist(:,loopi);%R位置处等效角度振幅
    Smove.f=Sinf.flist(:,loopi);%振动频率
    Smove.u=[zeros(1,Sinf.len);Smove.um*sin(2*pi*Smove.f*(0:Sinf.ts:((Sinf.len-1)*Sinf.ts)));zeros(1,Sinf.len)];
    Smove.R0=Sinf.Rlist(:,loopi);%理想时 主子间0时刻前初始相对位置
    Smove.Rf0=[0;0;0];%理想时 子0时刻前偏移
    Smove.Rf=[-Smove.R0(1)*sin(Smove.u(2,:)).*sin(Smove.u(2,:));zeros(1,Sinf.len);Smove.R0(1)*sin(Smove.u(2,:)).*cos(Smove.u(2,:))];
    Smove.R=Smove.R0+Smove.Rf;
    Smove.att=2*Smove.u;
    Smove.att_s0=[0;0;0]*arcdeg;%理想时 主子间0时刻前初始相对角度
    Smove.diffRf0=[0;0;0];%Rf的初始变化率
    Smove.U0=Smove.diffRf0+cross(Sinf.wim0,Smove.R0);
    Smove.atterr0=Sinf.aerrlist(:,loopi);%注意这里旋转顺序zxy
    SmoveCell{1,loopi}=Smove;
end