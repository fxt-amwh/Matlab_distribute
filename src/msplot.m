function msplot(mnp, x, y, xstr, ystr)
if mod(mnp,10)==1, figure; end   % ����ǵ�һ��Сͼ�����½�һ��figure
subplot(mnp); plot(x, y);  grid on;
if nargin==4, ystr = xstr; xstr = 't / s'; end  % ���ֻ����һ���ַ�������Ĭ��xlabelΪʱ��
xlabel(xstr); ylabel(ystr);