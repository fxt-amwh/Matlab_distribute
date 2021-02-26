function SmoveCell=my_nSmovePack(Sinf)
% �������ƣ�my_nSmovePack
% �������ܣ���������Ե����ӹ����������ݽ��
% ���룺Sinf      :�ӹ���Ϣ���þ���
% �����SmoveCell :�ӹ߽ṹ��Ԫ��
arcdeg = pi/180;
n=size(Sinf.Rlist,2);
SmoveCell=cell(1,n);
for loopi=1:n
    Smove.um=Sinf.ulist(:,loopi);%Rλ�ô���Ч�Ƕ����
    Smove.f=Sinf.flist(:,loopi);%��Ƶ��
    Smove.u=[zeros(1,Sinf.len);Smove.um*sin(2*pi*Smove.f*(0:Sinf.ts:((Sinf.len-1)*Sinf.ts)));zeros(1,Sinf.len)];
    Smove.R0=Sinf.Rlist(:,loopi);%����ʱ ���Ӽ�0ʱ��ǰ��ʼ���λ��
    Smove.Rf0=[0;0;0];%����ʱ ��0ʱ��ǰƫ��
    Smove.Rf=[-Smove.R0(1)*sin(Smove.u(2,:)).*sin(Smove.u(2,:));zeros(1,Sinf.len);Smove.R0(1)*sin(Smove.u(2,:)).*cos(Smove.u(2,:))];
    Smove.R=Smove.R0+Smove.Rf;
    Smove.att=2*Smove.u;
    Smove.att_s0=[0;0;0]*arcdeg;%����ʱ ���Ӽ�0ʱ��ǰ��ʼ��ԽǶ�
    Smove.diffRf0=[0;0;0];%Rf�ĳ�ʼ�仯��
    Smove.U0=Smove.diffRf0+cross(Sinf.wim0,Smove.R0);
    Smove.atterr0=Sinf.aerrlist(:,loopi);%ע��������ת˳��zxy
    SmoveCell{1,loopi}=Smove;
end