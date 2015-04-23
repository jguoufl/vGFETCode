%anti_dummy function used to convert Ne to quasi Fermi energy
%(Zhibin Ren 6-5-00)

function [y]=anti_dummy(x,fermi_flag,fermi_order)

if fermi_order==0
  if fermi_flag==0
    y=log(x);
  elseif fermi_flag==1
    y=log(exp(x)-1);
  end
elseif fermi_order==1/2
  if fermi_flag==0
    y=log(x);
  elseif fermi_flag==1
    y=log(x)+3.53553e-1*x-4.95009e-3*x.^2+1.48386e-4*x.^3-4.42563e-6*x.^4;
  end
end

