function phi = qq2phi(qpb, qnb)
qerr = qmul(qnb, qconj(qpb));
phi = q2rv(qerr);