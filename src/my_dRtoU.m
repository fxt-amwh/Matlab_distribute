function U=my_dRtoU(R,varargin)
% �������ƣ�my_dRtoU
% �������ܣ����������λ��ʸ��΢�ַ��̣���Ӧmy_UtodiffU
% ���룺R      :��ǰ����ϵ�����λ��
%       varargin:������΢�����������,��Ԫ�����          
%              ÿһ�����������δ��dR λ��΢��
%                                  wim ��������ڹ�������ϵ���ٶ�����������ϵ��ͶӰ
% �����dR:��ǰ����ϵ�����λ��΢��
dR=varargin{1};
wim=varargin{2};
U=dR+cross(wim,R);