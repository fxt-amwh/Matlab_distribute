function fs=my_dUtofs(U,varargin)
% �������ƣ�my_dUtofs
% �������ܣ���������Ե����ٶ�΢�ַ���,��my_UtodiffU��Ӧ
% ���룺U      :��ǰ����ϵ������ٶ�
%       varargin:������΢�����������,��Ԫ�����          
%              ÿһ�����������δ��Cms��������������任����Csm���֣�3*3
%                                  dU  �ٶ�΢��
%                                  fm  ���߱���
%                                  wim ��������ڹ�������ϵ�ٶ�����������ϵ��ͶӰ
% �����dU:��ǰ����ϵ������ٶ�΢��
Cms=varargin{1};
dU=varargin{2};
fm=varargin{3};
wim=varargin{4};
fs=Cms*(dU+fm+cross(wim,U));