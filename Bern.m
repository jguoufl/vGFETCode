function [y]=Bern(x) % evaluate Bernoull'function

y=zeros(length(x),1);
for ii=1:length(x)
    if abs(x(ii))>1e-5
        y(ii)=x(ii)/(exp(x(ii))-1);
    else
        y(ii)=1;
    end
end
