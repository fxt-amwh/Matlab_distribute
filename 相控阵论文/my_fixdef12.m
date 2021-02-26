function dz=my_fixdef12(xdir,xm,ym,Zmax,flag,varargin)
ddx=0;ddy=0;
if nargin==7
    ddx=varargin{1};
    ddy=varargin{2};
else
    
end
%注意暂时使用
global N;
ydir=(-floor(N/2):1:floor(N/2))*0.05;

[MList,NList]=meshgrid(xdir,ydir);

MList=MList+ddx;
NList=NList+ddy;
switch flag
    case 1
        dz=Zmax*(MList.^2)/(xm^2);
    case 2
        dz=Zmax*(MList.^2+NList.^2)/(xm^2+ym^2);
    case 3
        dz=Zmax*(MList.^2/(xm^2)-NList.^2/(ym^2));     
end