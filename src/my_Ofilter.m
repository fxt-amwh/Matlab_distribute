function Ok=my_Ofilter(Ok_1,I,varargin)
% �������ƣ�my_Ofilter
% �������ܣ���ͨ�˲�������
% ���룺Ok_1     :��һʱ���˲����
%       I        :��ǰʱ������
%       varargin{1} :��Ϊ������ʱ��ʾ a ������ʱ��ʾfc ����Ƶ��
%       varargin{2} :�������ts ��ϵ�������
% �����Ok:�˲����
if nargin==3
    a=varargin{1};
elseif  nargin==4
    fc=varargin{1};%�˲�������Ƶ��Hz
    ts=varargin{2};
    Tc=1/(2*pi*fc);
    a=ts/(Tc+ts);
else
    a=0.5;
end
Ok=(1-a)*Ok_1+a*I;