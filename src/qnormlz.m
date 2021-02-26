function qnb = qnormlz(qnb)
nm = qnb'*qnb;
if nm<1e-6,  qnb = [1; 0; 0; 0];  % 表示姿态的四元数，其模值应约为1
else        qnb = qnb/sqrt(nm);    end