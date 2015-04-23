%---------------------------------------------------------------------------
%       gmapja(x,y,z,zmin,zmax,V,Vlabs,xlabel,ylabel,title) 
%
%       GMAPJA makes a 2D pseudocolor/contour plot of the data z = f(x,y).
%
%     - x and y are two matrix arguments that specify the domain overwhich
%       the function is to be plotted where (x(i,j),y(i,j)) is the i,jth
%       point in the domain.  With two vector arguments replacing the two 
%       matrix arguments, length(x) = n and length(y) = m where 
%       [m,n] = size(Z).  In this case, the vertices of the surface patches 
%       are the triples (x(j), y(i), Z(i,j)).  Note that x corresponds to 
%       the columns of Z and y corresponds to the rows.
%     - z is the data to be plotted.  We have z(i,j) = f(x(i,j),y(i,j)) or
%       Z = f(x(j),y(i)).
%     - zmin is an optional argument specifying the minimum value of z to be 
%       plotted.
%     - zmax is an optional argument specifying the maximum value of z to be 
%       plotted.
%     - V is an optional argument which is either a positive integer or an 
%       array of values.  If it is a positive integer, it specifies the
%       number of contour lines to be drawn.  If it is an array, it specifies
%       that LENGTH(V) contour lines are to be drawn at the values specified 
%       in vector V.  To compute a single contour at the level v, 
%       use V = [v,v]  The default is to draw 10 automatically chosen contour 
%       lines.
%     - Vlabs is an optional argument which is either a positive integer or  
%       an array of values or the text string 'manual'.  A positive integer 
%       specifies the number of contour lines to be labelled with their values. 
%       An array specifies that LENGTH(Vlabs) contour lines of the values
%       specified in Vlabs are to be labelled with their values.  If a contour 
%       line specified in Vlabs does not appear on the plot, the function will 
%       attempt to find and label the nearest valued contour line.  The 
%       character string 'manual' allows placement of contour labels at the 
%       locations clicked on with a mouse.  Pressing the return key terminates
%       labeling.  The default is to label all contour lines.
%     - xlabel is an optional x-axis label.  Default is blank.
%     - ylabel is an optional y-axis label.  Default is blank.
%     - title is an optional title.  Default is blank.

%-----------------------------------------------------------------------------

function gmapja(X,Y,Z,zmin1,zmax1,V1,Vlabs1,xlabel1,ylabel1,title1)

%
%  verify input arguments
%
error(nargchk(3,10,nargin));
x = X;
y = Y;
z = Z;
k = find(z~=-Inf & z~=Inf & z~=NaN);
zchop = z(k);
if nargin < 4,
   zmin = min(min(zchop));
else
   zmin = zmin1;
   if isempty(zmin), zmin = min(min(zchop)); end
end
if nargin < 5,
   zmax = max(max(zchop));
else
   zmax = zmax1;
   if isempty(zmax), zmax = max(max(zchop)); end
end
if nargin < 6,
   ncntrs = 10;
   dz = (zmax-zmin)/(ncntrs+1);
   V = zeros(1,ncntrs);
   for icntr=1:ncntrs,
       V(icntr) = zmin + icntr*dz;
   end
else
   if prod(size(V1)) == 0,
      ncntrs = 10;
      dz = (zmax-zmin)/(ncntrs+1);
      V = zeros(1,ncntrs);
      for icntr=1:ncntrs,
          V(icntr) = zmin + icntr*dz;
      end
   elseif prod(size(V1)) == 1,
      ncntrs = V1;
      dz = (zmax-zmin)/(ncntrs+1);
      V = zeros(1,ncntrs);
      for icntr=1:ncntrs,
          V(icntr) = zmin + icntr*dz;
      end
   else
      V = V1;
   end      
end
if nargin < 7,
   Vlabs = V;
else
   if prod(size(Vlabs1)) == 0,
      Vlabs = V;
   elseif prod(size(Vlabs1)) == 1,
      nlabels = Vlabs1;
      dzlab = (zmax-zmin)/(nlabels+1);
      Vlabs = zeros(1,nlabels);
      for ilabel=1:nlabels,
          Vlabs(ilabel) = zmin + ilabel*dzlab;
      end
   else
      Vlabs = Vlabs1;
   end      
end
if nargin < 8,
   xlab = '';
else
   xlab = xlabel1;    
   if ~isstr(xlab),
      error('xlabel must be a string.')
   end
end
if nargin < 9,
   ylab = '';
else
   ylab = ylabel1;    
   if ~isstr(ylab),
      error('ylabel must be a string.')
   end
end
if nargin < 10,
   tit = '';
else
   tit = title1;     
   if ~isstr(tit),
      error('title must be a string.')
   end
end

%
%  chop the data
%
zchop = z;
k1 = find(zchop<zmin);
zchop(k1) = zmin;
k2 = find(zchop>zmax);
zchop(k2) = zmax;

%
%  draw the pseudocolor plot
%
colormap('jet')
pcolor(x,y,zchop)
caxis([zmin zmax])
shading('interp')
hold on

%
%  draw the contour plot
%
   if prod(size(V))==1,
      V = [V,V];
   end
   [cs,h] = contour(x,y,zchop,V,'k');

%
%  label appropriate contours
%
   nlabels = length(Vlabs);
   for ilabel=1:nlabels,
      [errmin,kcntr] = min(abs(V-Vlabs(ilabel)));
      Vlabs(ilabel) = V(kcntr);
  end
  clabel(cs);
  
   hold off

%
%  add colorbar
%
h1=colorbar('v');
set(h1, 'fontsize',[28]);
%
%  label axes
%
%set(gca,'YDir','reverse');
xlabel(xlab)
ylabel(ylab)
title(tit)
set(gca, 'fontsize',[28], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.65 0.70]);
