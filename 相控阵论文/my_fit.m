function [xi,yi,varargout]=my_fit(xinput,yinput,varargin)
%功能：拟合增益衰减曲线
dN=1000;%细分度
if nargin==3
    p0=varargin{1};
else
%     p0=[0,0,1];
    p0=[0,1];
% p0=[1];
end
sym t;
fuc=fittype('t*exp(a)*(exp(k*t)-1)','independent','t','coefficients',{'a','k'});
% fuc=fittype('t*(exp(k*t)-1)','independent','t','coefficients',{'k'});
cfun=fit(xinput',yinput',fuc,'StartPoint',p0);
xi=linspace(min(xinput),max(xinput),dN);
yi=cfun(xi);
if nargout==3
    varargout{1}=cfun;
end