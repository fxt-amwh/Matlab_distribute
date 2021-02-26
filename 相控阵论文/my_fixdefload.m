function varargout=my_fixdefload(xdir,L,u,varargin)
flag=0;
dz=0;
if 4==nargin
    flag=varargin{1};
end
if 0==flag
    dz=my_xtodz(xdir,L,u);
elseif 1==flag %Í¨¹ı»úÒíÒí¼âÎ»ÖÃ£¬Òí¼âÄÓÇú½Ç
    [dx,dz]=my_XstdToDzMap(xdir,L,u);
elseif 2==flag
    [dx,dz]=my_XstdToDzMap(xdir,L,u,1e-5,flag);
else
end
if 1==nargout
    varargout{1}=dz;
end
if 2==nargout
    varargout{1}=dx;
    varargout{2}=dz;
end