%% һ��+���ӣ��߶�ģ�黯��
% 1��������ƹ켣�˶����������ߴ��������
% 2��������ƻ����˶��������ӹߴ���������
% 3���ӹ�����ע�����
% 4����������Ե���
% 5���˲�
% 6��������
%% �����˶�������Ʒ���
clc
clear
close all
gvar;    % ����ȫ�ֱ���
global setPara;
d2g=180/pi;
g2d=1/d2g;
setPara.ts=0.01;%����ʱ��
setPara.att0=[0;0;0]*d2g;%�ػ���ʼ��̬
setPara.vn0=[0;200;0];%�ػ���ʼ�ٶ�
setPara.latitude=34*d2g;%γ��
setPara.longitude=108*d2g;%����
setPara.h=105;%��
              %����������   ��ת������   ���������  ���ٶ�     ʱ��
              %  ��/s      ��/s        ��/s        m/s2       s
setPara.wat=[    0,         0,          0,         0,         10       %��ֹ
                 1,         0,          0,         0,         60       %����
                -1,         0,          0,         0,         60       %����
                 0,         0,          0,         0,         50       %����
                 0,         0,          1,         0,         360      %����ת��
                 0,         0,          0,         0,         50       %����
                -1,         0,          0,         0,         60       %����
                 1,         0,          0,         0,         60       %����
                 0,         0,          0,         0,         10       %����
        ];    %��ֹ

ts = setPara.ts;%����ʱ��
att0 = setPara.att0; vn0 = setPara.vn0; pos0 = [[setPara.latitude;setPara.longitude];setPara.h];
wat = setPara.wat;
wat(:,1:3) = wat(:,1:3)*arcdeg;  % deg/min->deg/s
fprintf('�켣������...'); 
[attm, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
fprintf('�켣������ɣ�\n���߷�����...');
[wm, vm] = my_av2imu(attm, vn, pos, ts);
tt = (0:length(attm)-1)'*ts;
% �켣��ͼ
msplot(221, tt, attm/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(222, tt, vn, 'Vel / m.s^{-1}'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(223, tt, deltapos(pos), '\DeltaPos / m');
      legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')
msplot(224, pos(:,2)/arcdeg, pos(:,1)/arcdeg, '\itL\rm / ( \circ )', '\it\lambda\rm / ( \circ)');
      hold on, plot(pos(1,2)/arcdeg, pos(1,1)/arcdeg, 'ro');
%  ���Դ�������Ϣ��ͼ
msplot(121, tt(2:end), wm/ts/arcdeg, '\it\omega^b_{ib}\rm / ( \circ.s^{-1} )');
      legend('\it\omega^b_{ibx}', '\it\omega^b_{iby}', '\it\omega^b_{ibz}');
msplot(122, tt(2:end), vm/ts, '\itf^b\rm_{sf} / ( m.s^{-2} )');
      legend('\itf^b\rm_{sf\itx}', '\itf^b\rm_{sf\ity}', '\itf^b\rm_{sf\itz}');
% ��ά�켣
figure
delpos=deltapos(pos);
plot3(delpos(:,1),delpos(:,2),delpos(:,3));
% axis([0,-8000,]);
grid on;
xlabel('γ�Ⱦ���/m'),ylabel('���Ⱦ���/m'),zlabel('��/m');
%% ������Ϣ
wim0=my_getII(wm(1,:)',ts);
wim=[wim0,my_getII(wm',ts)];%wimע������������ϵ�µ�����
fm=[my_getII(vm(1,:)',ts),my_getII(vm',ts)];
MINS.wim=wim;MINS.wim0=wim0;MINS.fm=fm;MINS.ts=ts;
fprintf('���߷�����ɣ�\n');