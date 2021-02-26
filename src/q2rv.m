function rv = q2rv(q) 
if q(1)<0,  q = -q;  end
nmhalf = acos(q(1));  % 等效旋转矢量模值的一半
if nmhalf>1e-20,  b = 2*nmhalf/sin(nmhalf);
else            b = 2;                end
rv = b*q(2:4);