function dvul=my_relgetdvul(Csm,Vm,Rm,wm_m, vm_m,ws_m,vs_m,ts)
nts = size(wm_m,1)*ts;%单次更新时间间隔
a=(3*(wm_m(1,:))-(wm_m(2,:)));
b=(wm_m(1,:))-(wm_m(2,:));
wibm_m=a+2*b*nts;
dvul=-(cross(wibm_m,Vm)+2*cross(b,Rm)+cross(wibm_m,cross(wibm_m,Rm)))*nts;