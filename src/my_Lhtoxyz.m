function xyz=my_Lhtoxyz(Lh)
L=Lh(1);
a=Lh(2);
h=Lh(3);
Re = 6.378136998405e6;
ff = 1/298.257223563; ee = sqrt(2*ff-ff^2); e2 = ee^2; 
RN=Re/sqrt(1-e2*(sin(L)^2));
x=(RN+h)*cos(L)*cos(a);
y=(RN+h)*cos(L)*sin(a);
z=(RN*(1-e2)+h)*sin(L);
xyz=[x,y,z];
