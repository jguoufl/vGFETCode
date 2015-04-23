%fermi means integral, the integral order is indicated by the 
%last input argument
%x is the input varible, and fermi_flag is the degeneracy flag 
%1 for degenerate calculation, 0 for nondegenerate calculation

 
function [y]=fermi(x,fermi_flag,fermi_order)

if fermi_order==1/2
   
if fermi_flag==1
  exp_fac=exp(-0.17*(x+1.0).^2);
  nu=x.^4+50.0+33.6*x.*(1.0-0.68*exp_fac);
  zeta=3.0*sqrt(pi)./(4.0*nu.^0.375);
  y=exp(x)./(1.0+zeta.*exp(x));	
elseif fermi_flag==0
  y=exp(x);
end   

elseif fermi_order==0

if fermi_flag==1
  y=log(1+exp(x));
elseif fermi_flag==0
  y=exp(x);
end

elseif fermi_order==-1/2
   
delta_x=1e-8;
if fermi_flag==1
   x_r=x+delta_x;
   exp_fac=exp(-0.17*(x_r+1.0).^2);
   nu=x_r.^4+50.0+33.6*x_r.*(1.0-0.68*exp_fac);
   zeta=3.0*sqrt(pi)./(4.0*nu.^0.375);
   y_r=exp(x_r)./(1.0+zeta.*exp(x_r));	
   
   x_l=x-delta_x;
   exp_fac=exp(-0.17*(x_l+1.0).^2);
   nu=x_l.^4+50.0+33.6*x_l.*(1.0-0.68*exp_fac);
   zeta=3.0*sqrt(pi)./(4.0*nu.^0.375);
   y_l=exp(x_l)./(1.0+zeta.*exp(x_l));	

   y=(y_r-y_l)./(2*delta_x);
elseif fermi_flag==0
   y=exp(x);
end
   
end   
