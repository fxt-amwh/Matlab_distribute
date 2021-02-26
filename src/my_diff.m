function dx=my_diff(x,varargin)
% 函数名称：my_diff
% 函数功能：反解算输入序列微分(一阶微分)
% 输入：x      :输入序列
%       f      :功能选择 
%               微分方向 1 行微分 2列微分
%               微分初值
%               微分补偿模式 1 不补偿 2 线性补偿
% 输出：dx:微分
switch nargin
    case 1
        if size(x,1)>1
            dir=1;
        else
            dir=2;
        end
        flag=1;
    case 2
        dir=varargin{1};
        flag=1;
    case 3
        dir=varargin{1};
        x0=varargin{2};
        flag=1;
        if dir==1
            x=[x0;x];
        else
            x=[x0,x];
        end
        
    case 4
        dir=varargin{1};
        x0=varargin{2};
        flag=varargin{3};
        if dir==1
            x=[x0;x];
        else
            x=[x0,x];
        end
end

if flag==1
    dx=diff(x,1,dir);
elseif flag==2
    dx=diff(x,1,dir);
    if length(dx)<=1
        
    else
        ddx=diff(dx,1,dir);
        if dir==1
            dx(1:end-1,:)=dx(1:end-1,:)+0.5*ddx;
            dx(end,:)=dx(end,:)+0.5*ddx(end,:);
        elseif dir==2
            dx(:,1:end-1)=dx(:,1:end-1)+0.5*ddx;
            dx(:,end)=dx(:,end)+0.5*ddx(:,end);
        end
    end
end

