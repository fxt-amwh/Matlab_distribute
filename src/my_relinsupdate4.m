function [qms,vn,U,R] = my_relinsupdate4(qms,U,R,wm_m, vm_m, ws_m, vs_m, ts)
wm_arg=my_getII(wm_m',ts);
fm_arg=my_getII(vm_m',ts);
ws_arg=my_getII(ws_m',ts);
fs_arg=my_getII(vs_m',ts);
wm_3=[wm_arg(:,1)*3/2-wm_arg(:,2)/2,wm_arg(:,1)/2+wm_arg(:,2)/2,...
    wm_arg(:,2)*3/2-wm_arg(:,1)/2];
fm_3=[fm_arg(:,1)*3/2-fm_arg(:,2)/2,fm_arg(:,1)/2+fm_arg(:,2)/2,...
    fm_arg(:,2)*3/2-fm_arg(:,1)/2];
ws_3=[ws_arg(:,1)*3/2-ws_arg(:,2)/2,ws_arg(:,1)/2+ws_arg(:,2)/2,...
    ws_arg(:,2)*3/2-ws_arg(:,1)/2];
fs_3=[fs_arg(:,1)*3/2-fs_arg(:,2)/2,fs_arg(:,1)/2+fs_arg(:,2)/2,...
    fs_arg(:,2)*3/2-fs_arg(:,1)/2];
phim = my_getdphim(wm_m);
phis = my_getdphim(ws_m);
qms_1=qms;
qms_2 = qmul(rv2q(-phis/2), qmul(qms, rv2q(phim/2)));
qms = qmul(rv2q(-phis), qmul(qms, rv2q(phim)));  % 姿态更新
qms = qnormlz(qms);
qms_3=[qms_1,qms_2,qms];
u{1}=cell(4,1);
u{2}=cell(4,1);
u{3}=cell(4,1);
for i=1:3
    u{i}{1}=q2mat(qms_3(:,i));
    u{i}{2}=fs_3(:,i);
    u{i}{3}=fm_3(:,i);
    u{i}{4}=wm_3(:,i);
end
U1=U;
U=my_LK(@my_UtodiffU,U,u,2*ts,4);% 速度更新
U_3=[U1,(U+U1)/2,U];
u1{1}=cell(2,1);
u1{2}=cell(2,1);
u1{3}=cell(2,1);
for i=1:3
    u1{i}{1}=U_3(:,i);
    u1{i}{2}=wm_3(:,i);
end
vn=U-cross(wm_3(:,3),R);
R=my_LK(@my_RtodiffR,R,u1,2*ts,4);% 位置更新