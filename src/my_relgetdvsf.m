function dvsf=my_relgetdvsf(Csm,wm_m, vm_m,ws_m,vs_m)
dvsfm_m=sum(vm_m,1);
dvsfs_m=sum(vs_m,1);
dwm_m=sum(wm_m,1);
dws_m=sum(ws_m,1);
dvsf_2=0.5*cross(dwm_m,dvsfs_m*(Csm)')+...
    2/3*(cross(wm_m(1,:),vs_m(2,:)*(Csm)')+...
    cross(vs_m(1,:)*(Csm)',wm_m(2,:)));
dvsf_3=0.5*cross(dws_m,dvsfs_m)+...
    2/3*(cross(ws_m(1,:),vs_m(2,:))+cross(vs_m(1,:),ws_m(2,:)));
dvsf=-dvsfm_m+dvsf_2+(dvsfs_m+dvsf_3)*(Csm)';