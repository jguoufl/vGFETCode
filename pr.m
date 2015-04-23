% generate 3D figures
function [xlin, ylin, Z]=pr(n)

x=n(:,1)/1e-9;y=n(:,2)/1e-9;z=n(:,3);   % convert to nm
nombx=500  
nomby=500  
xlin = linspace(min(x),max(x),nombx);
ylin = linspace(min(y),max(y),nomby);
[X,Y] = meshgrid(xlin,ylin);
Z = griddata(x,y,z,X,Y);
viewmode=2;
if viewmode==1
    surf(X,Y,Z)
    shading interp;
    axis tight;
    %view(0,90)
    colorbar
    colormap(0.9*jet+0.1*flag)
else
    gmapja(X,Y,Z);
    h_xlabel=get(gca, 'xlabel');    h_ylabel=get(gca, 'ylabel');
    set(h_xlabel,'string','x [nm]','fontsize',[28]);
    set(h_ylabel,'string','y [nm]','fontsize',[28]);
end
xlin=xlin*1e-9;     % convert to m
ylin=ylin*1e-9;     % convert to m


