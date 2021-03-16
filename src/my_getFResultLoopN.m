function [INSSFR,Filter,NINSS]=my_getFResultLoopN(MINS,SINS,KFinit,atterr,len,selcCell,flag)
% 函数名称：my_getFResultLoopN 
% 函数功能：由子惯输出进行相对导航滤波
% 输入：MINS      :主惯
%       SINS      :子惯
%       KFinit    :子惯滤波器初始值
%       atterr    :子惯安装误差角估计初始值
%       len       :滤波长度       
%       selcCell  :select=0，selcCell无用 接口获得无误差数据
%                  select=1，selcCell表示 SINSR 子惯纯相对导航结果
%                  select=2，selcCell表示 SERR  子惯误差初始值
%       select    :选择子惯误差来源
% 输出：INSSFR:每个子惯结果元组
%       Filter:每个子惯滤波器中间值元组
%       NINSS :融合后最终结果结构体
fprintf('滤波...%5.0f %%',0);
[M,N]=size(SINS);%预留升级面阵SINS接口 M=1为线SINS组
gvar;    % 加载全局变量
looplen=floor(len);
ts=SINS{1,1}.ts;%采样时间
nts=2*ts;%解算时间
%初始分配
Rloop=cell(M,N);%相对导航 R 递推变量
Uloop=cell(M,N);%相对导航 U 递推变量
INSSFR=cell(M,N);%相对导航滤波结果
Filter=cell(M,N);%相对导航滤波中间值
vn=cell(M,N);%相对导航滤波速度中间值
qms=cell(M,N);%相对导航 q 递推变量
attnow=cell(M,N);%滤波 角度 递推变量
dR=cell(M,N);%滤波量测 dR 递推变量
attint=cell(M,N);%初始角 用来计算角度误差 旋转次序yxz
kfft=cell(M,N);%滤波器配套状态转移矩阵
kf=cell(M,N);%卡尔曼滤波器
RFilter=cell(M,N);%低通滤波器初始值
ws_m=cell(M,N);%循环时每个子惯的两次量测 角速度
fs_m=cell(M,N);%循环时每个子惯的两次量测 比力

NINSS.uf=zeros(3,looplen);
for m=1:M
    for n=1:N
        Rloop{m,n}=SINS{m,n}.R0;
        Uloop{m,n}=cross(MINS.wim0,SINS{m,n}.R0);

        INSSFR{m,n}.Rall=zeros(looplen,3);
        INSSFR{m,n}.vnall=zeros(looplen,3);
        INSSFR{m,n}.qmsall=zeros(looplen,4);
        INSSFR{m,n}.attall=zeros(3,looplen);

        Filter{m,n}.X=zeros(15,looplen);
        
        Cmserr_1=my_a2mat(atterr{m,n});%注意中间处理旋转顺序yxz 不准确的初始角
        qms{m,n}=m2qua(Cmserr_1);%初始姿态 不准确的初始角
        attint{m,n}=m2att(Cmserr_1);%旋转次序yxz 输出依次俯仰角 横滚角 方位角

        fs0=[0;0;0];
        kfft{m,n}=my_kfft15(MINS.wim0,q2mat(qms{m,n})',fs0,nts);
        kf{m,n}= kfinit(KFinit{m,n}.Qk, KFinit{m,n}.Rk, KFinit{m,n}.P0,kfft{m,n}.phi,kfft{m,n}.H);  % kf滤波器初始化
        kf{m,n}.Xk(10:12) = [0.1;0.1;0.1]*dph;
        kf{m,n}.Xk(13:15) = [20;20;20]*ug;
        RFilter{m,n}=SINS{m,n}.R0;
    end
end

for i=1:looplen
    if flag.select==0%无误差
        for m=1:M
            for n=1:N
                ws_m{m,n}=SINS{m,n}.wis(:,(2*i-1):(2*i))';
                fs_m{m,n}=SINS{m,n}.fs(:,(2*i-1):(2*i))';
            end
        end
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; 
        fm_m=MINS.fm(:,(2*i-1):(2*i))';
    elseif flag.select==1%直接使用相对导航的误差
        SINSR=selcCell;
        for m=1:M
            for n=1:N
            ws_m{m,n}=SINSR{m,n}.ws_m_addnoise(:,:,i);
            fs_m{m,n}=SINSR{m,n}.fs_m_addnoise(:,:,i);
            end
        end
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; 
        fm_m=MINS.fm(:,(2*i-1):(2*i))';
    else %自注入误差
        SERR=selcCell;
        for m=1:M
            for n=1:N
                [ws_m{m,n}, fs_m{m,n}] = my_imuadderr...
                    (SINS{m,n}.wis(:,(2*i-1):(2*i))',...
                    SINS{m,n}.fs(:,(2*i-1):(2*i))', ...
                    SERR{m,n}.eb,SERR{m,n}.web,...
                    SERR{m,n}.db,SERR{m,n}.wdb,ts);%子惯注入噪声
            end
        end
        wm_m=MINS.wim(:,(2*i-1):(2*i))'; 
        fm_m=MINS.fm(:,(2*i-1):(2*i))';
    end
%—————————————————————————————————————————————————
% 以下可以获得每次时刻的主惯（1个）和子惯（M*N个）传感器量测
% 量测为 wm_m 两次主惯角速率 fm_m 两次主惯比力
%        ws_m 两次子惯角速率 fs_m 两次子惯比力 以元组{M,N}使用
% 这里为了后续多节点关联方法设计使用
    for m=1:M
        for n=1:N
            [qms{m,n},vn{m,n},Uloop{m,n},Rloop{m,n}] = my_relinsupdate5...
                (qms{m,n},Uloop{m,n},Rloop{m,n},wm_m,fm_m,ws_m{m,n},fs_m{m,n},ts);% 子惯纯相对导航
            RFilter{m,n}=my_Ofilter(RFilter{m,n},Rloop{m,n},2,nts);%低通滤波
        end
    end
    if flag.EnFusion  
        for m=1:M
            for n=1:N
                attnow{m,n}=q2att(qms{m,n});
                uf=attnow{m,n}-attint{m,n};%误差角

                dR{m,n}=my_getdR(RFilter{m,n},SINS{m,n}.R0,uf);%求出量测

                wim=sum(wm_m)/2;%可能有问题
                fs=sum(fm_m)/2;
                kfft{m,n}=my_kfft15(wim',q2mat(qms{m,n})',fs',nts);
                kf{m,n}.Phikk_1=kfft{m,n}.phi;
                kf{m,n}.Gammak=kfft{m,n}.Gammak;

                kf{m,n} = kfupdate(kf{m,n},dR{m,n},'B');%卡尔曼滤波

                Filter{m,n}.X(:,i)=kf{m,n}.Xk;
                
                attint{m,n}=q2att(qdelphi(a2qua(attint{m,n}),-kf{m,n}.Xk(1:3)));% 2021 03 15新增 更新初始姿态误差阵
                
                qms{m,n} = qdelphi(qms{m,n},-kf{m,n}.Xk(1:3));
                kf{m,n}.Xk(1:3) = 0;  % ·反馈
                Uloop{m,n}= Uloop{m,n}- kf{m,n}.Xk(4:6);  kf{m,n}.Xk(4:6) = 0;
                Rloop{m,n}= Rloop{m,n}- kf{m,n}.Xk(7:9);  kf{m,n}.Xk(7:9) = 0;

                INSSFR{m,n}.qmsall(i,:)=qms{m,n}';
                INSSFR{m,n}.Rall(i,:)=Rloop{m,n}';
                INSSFR{m,n}.vnall(i,:)=vn{m,n}';

                [rz,rx,ry]=dcm2angle(q2mat(INSSFR{m,n}.qmsall(i,:)),'zxy');%转换到传统的坐标
                INSSFR{m,n}.attall(:,i)=[rx,ry,-rz];
            end
        end
        qms_all=my_CellAvr(qms);
        Uloop_all=my_CellAvr(Uloop);
        vn_all=my_CellAvr(vn);
        Rloop_all=my_CellAvr(Rloop);
        Rfall=Rloop_all-SINS{m,n}.R0;
        if Rfall(3)>0
            NINSS.uf(:,i)=[0;asin(sqrt(abs(Rfall'*Rfall)/abs(SINS{m,n}.R0'*SINS{m,n}.R0)));0];
        else
            NINSS.uf(:,i)=[0;-asin(sqrt(abs(Rfall'*Rfall)/abs(SINS{m,n}.R0'*SINS{m,n}.R0)));0];
        end
        
        for m=1:M % 最小二乘反馈 没有更新方差阵
            for n=1:N
%                 kf{m,n}.Xk(1:3) = q2rv(qdelphi(qms{m,n},qms_all));
%                 qms{m,n}= qms_all;  
%                 kf{m,n}.Xk(4:6) = Uloop{m,n}-Uloop_all;
%                 Uloop{m,n}= Uloop_all;
                kf{m,n}.Xk(7:9) = Rloop{m,n}-Rloop_all;
                Rloop{m,n}= Rloop_all;
            end
        end
        
        NINSS.qmsall(i,:)=qms_all';
        NINSS.Rall(i,:)=Rloop_all';
        NINSS.vnall(i,:)=vn_all';
        attMat=a2mat(2*NINSS.uf(:,i));
        [rz,rx,ry]=dcm2angle(attMat,'zxy');%转换到传统的坐标
        NINSS.attnewall(:,i)=[rx,ry,-rz];
        [rz,rx,ry]=dcm2angle(q2mat(NINSS.qmsall(i,:)),'zxy');%转换到传统的坐标
        NINSS.attall(:,i)=[rx,ry,-rz];
        
    else
        
        for m=1:M
            for n=1:N
                attnow{m,n}=q2att(qms{m,n});
                uf=attnow{m,n}-attint{m,n};%误差角

                dR{m,n}=my_getdR(RFilter{m,n},SINS{m,n}.R0,uf);%求出量测

                wim=sum(wm_m)/2;%可能有问题
                fs=sum(fm_m)/2;
                kfft{m,n}=my_kfft15(wim',q2mat(qms{m,n})',fs',nts);
                kf{m,n}.Phikk_1=kfft{m,n}.phi;
                kf{m,n}.Gammak=kfft{m,n}.Gammak;

                kf{m,n} = kfupdate(kf{m,n},dR{m,n},'B');%卡尔曼滤波

                Filter{m,n}.X(:,i)=kf{m,n}.Xk;

            %     [rz,rx,ry]=dcm2angle(q2mat(rv2q(kf{m,n}.Xk(1:3))),'zxy');
            %     Filter{m,n}.atterr=[rx;ry;-rz];

            %     kf{m,n}.Xk(1:3) = 0;
                qms{m,n} = qdelphi(qms{m,n},-kf{m,n}.Xk(1:3));
                kf{m,n}.Xk(1:3) = 0;  % ·反馈
                Uloop{m,n}= Uloop{m,n}- kf{m,n}.Xk(4:6);  kf{m,n}.Xk(4:6) = 0;
                Rloop{m,n}= Rloop{m,n}- kf{m,n}.Xk(7:9);  kf{m,n}.Xk(7:9) = 0;

                INSSFR{m,n}.qmsall(i,:)=qms{m,n}';
                INSSFR{m,n}.Rall(i,:)=Rloop{m,n}';
                INSSFR{m,n}.vnall(i,:)=vn{m,n}';

                [rz,rx,ry]=dcm2angle(q2mat(INSSFR{m,n}.qmsall(i,:)),'zxy');%转换到传统的坐标
                INSSFR{m,n}.attall(:,i)=[rx,ry,-rz];
            end
        end
        
        NINSS=INSSFR{1,1};%最终融合返回值 暂时不融合
    end
%————————————————————————————————————————————————————————
    if mod(i,looplen/100)==0
              fprintf('\b\b\b\b\b\b\b%5.0f %%', i/looplen*100);			%进度显示
    end
end
fprintf('\n滤波完成！\n'); 