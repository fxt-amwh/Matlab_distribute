function dv=my_relgetdv(Csm,Vm,Rm,wm_m, vm_m,ws_m,vs_m,ts)
dvsf=my_relgetdvsf(Csm,wm_m, vm_m,ws_m,vs_m);
dvul=my_relgetdvul(Csm,Vm,Rm,wm_m, vm_m,ws_m,vs_m,ts);
dv=dvsf+dvul;
