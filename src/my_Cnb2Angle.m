function angle=my_Cnb2Angle(Cnb)
pitch = asin(Cnb(3,2));
if abs(Cnb(3,2))<=0.999999
   roll = -atan2(Cnb(3,1),Cnb(3,3));
   yaw = atan2(Cnb(1,2),Cnb(2,2));
else
   roll = -atan2(Cnb(3,1),Cnb(3,3));
   yaw = 0;
end
T22=Cnb(2,2);
T33=Cnb(3,3);
T12=Cnb(1,2);
if abs(T22)<0.000001
    if T12>0
        yaw=pi/2;
    else
        yaw=-pi/2;
    end
elseif T22<0
    if T12>0
        yaw=yaw+pi;
    else
        yaw=yaw-pi;
    end
end
if roll>0&&T33<0
    roll=roll-pi;
elseif T33<0
    roll=roll+pi;
end
angle=[pitch,roll,yaw];