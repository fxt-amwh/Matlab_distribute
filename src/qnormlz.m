function qnb = qnormlz(qnb)
nm = qnb'*qnb;
if nm<1e-6,  qnb = [1; 0; 0; 0];  % ��ʾ��̬����Ԫ������ģֵӦԼΪ1
else        qnb = qnb/sqrt(nm);    end