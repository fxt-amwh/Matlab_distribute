function Cnb = my_a2mat(att)
s = sin(att); c = cos(att);
s1 = s(1); s2 = s(2); s3 = s(3);   c1 = c(1); c2 = c(2); c3 = c(3);
Cnb = [c2*c3+s2*s3*s1 -c2*s3+s2*c3*s1 -s2*c1
    s3*c1 c3*c1 s1
    s2*c3-c2*s3*s1 -s2*s3-c2*c3*s1 c2*c1];