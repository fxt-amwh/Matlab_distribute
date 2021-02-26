function [wm, fm] = my_imuadderr(wm, fm, eb, web, db, wdb, ts)
%wm:
%vm:
%eb:陀螺零偏arcdeg/s
%web:角度随机游走arcdeg/sqrt(s)
%db:
%wdb:
%ts:
m = size(wm,1); sts = sqrt(ts);
wm = wm + [ eb(1) + sts*web(1)*randn(m,1)/ts,eb(2) + sts*web(2)*randn(m,1)/ts,eb(3) + sts*web(3)*randn(m,1)/ts ];
fm = fm + [ db(1) + sts*wdb(1)*randn(m,1)/ts,db(2) + sts*wdb(2)*randn(m,1)/ts,db(3) + sts*wdb(3)*randn(m,1)/ts ];