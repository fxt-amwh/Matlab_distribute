function Ok=my_Ofilter(Ok_1,I,varargin)
% 函数名称：my_Ofilter
% 函数功能：低通滤波递推器
% 输入：Ok_1     :上一时刻滤波输出
%       I        :当前时刻输入
%       varargin{1} :当为三输入时表示 a 四输入时表示fc 带宽频率
%       varargin{2} :采样间隔ts 配合第三输入
% 输出：Ok:滤波输出
if nargin==3
    a=varargin{1};
elseif  nargin==4
    fc=varargin{1};%滤波器带宽频率Hz
    ts=varargin{2};
    Tc=1/(2*pi*fc);
    a=ts/(Tc+ts);
else
    a=0.5;
end
Ok=(1-a)*Ok_1+a*I;