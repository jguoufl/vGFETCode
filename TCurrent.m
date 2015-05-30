%%% Fuunction for Calculating Tunenling Current%%%%%%%%%%%
%%%%%%%%%%%%%%%Wenchao Chen, Sept,2012%%%%%%%%%%%%%%%%%%%

function [J_total]=TCurrent(Ec3D,Fn3D,Gnode1,NodeA,nox,N2D,phid,z_gc,nloc,Nz)

%%%z_gc: graphene+channel, grid in z direction
global m0 
global hbar q 

T=300;
kB=1.38e-23;
Effectivemass=0.5;
Mass=Effectivemass*m0;
Richardson=4*pi*q*Effectivemass*m0*kB^2/(2*pi*hbar)^3;

z=z_gc(1:nloc+1);                             %%% local mesh grid

for ii_n=1:N2D
    
    if isempty(find(Gnode1==ii_n))==1
        Jt(ii_n)=0;
        Jte(ii_n)=0;
    else
        Ec=Ec3D(nox+1:nloc+nox+1,ii_n);         %%%band profile in local mesh grid
        Fn=Fn3D(nox+1:Nz,ii_n);                 
        Fn_end=Fn(length(Fn));                  %%%quasi fermi level at drain side
        phis=Ec(1);                             %%%shottky barrier height
        for ii_z=1:length(z)
            if ii_z==1
                Tr(ii_z)=1;
            else
                Tr(ii_z)=exp(-1*abs(2/hbar*trapz(z(1:ii_z),abs(sqrt(2*Mass*q*(phis-Ec(1:ii_z)))))));
            end
        end
        dE=Ec(1:length(Ec)-1)-Ec(2:length(Ec));
        T_bridge=q*Richardson*T/kB*(log((1+exp((-Ec)*q/kB/T)))-log(1+exp((Fn_end-Ec)*q/kB/T)));%%%For tunneling
        Jt(ii_n)=abs(sum(Tr(1:length(Tr)-1)'.*dE.*T_bridge(1:length(T_bridge)-1)));       %%tunneling current
        
        N_grid=200;
        EnergyGrid=linspace(0,0.5,N_grid);
        dEnergy=EnergyGrid(2)-EnergyGrid(1);
        TE_bridge=q*Richardson*T/kB*(log((1+exp((-phis-EnergyGrid)*q/kB/T)))-log(1+exp((Fn_end-phis-EnergyGrid)*q/kB/T)));
        Jte(ii_n)=abs(dEnergy*sum(TE_bridge));                           %%thermionic current
    end
   
end
% J_total=sum((Jt+Jte)'.*NodeA)/sum(NodeA);
% J_tunnel=sum((Jt)'.*NodeA)/sum(NodeA);
J_total=sum((Jt)'.*NodeA)/sum(NodeA);
    
    