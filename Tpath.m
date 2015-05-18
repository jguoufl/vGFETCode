%%% Fuunction for Searching for the tunneling path%%%%%%%
%%%%%%%%%%%%%%%Wenchao Chen, Sept,2012%%%%%%%%%%%%%%%%%%%

function [J_total,Jt]=Tpath(Ec3D,Fn3D,xlEcz, ylEcz, Ecz2D,Nodex,Nodey,Gnode1,NodeA,nox,N2D,phid,z_gc,nloc,Nz,Vg_bias)

%%% output: Jt, arry in the size of 2D nodes, tunneling current density at the nodes
%%% J_total, the tunneling current density, a number

global m0 
global hbar q 

L_side=max(Nodex);

T=300;
kB=1.38e-23;
Effectivemass=0.5;
Mass=Effectivemass*m0;
Richardson=4*pi*q*Effectivemass*m0*kB^2/(2*pi*hbar)^3;

z=z_gc(1:nloc+1);                             %%% local mesh grid

[X,Y] = meshgrid(xlEcz,ylEcz);

if Vg_bias>0   % the selection of the theta range for finding the tunnelng path, emperical
    theta=linspace(0,pi/2,30);                    %%%for choosing nonlocal mesh
else
    theta=linspace(pi-pi/40,pi/2,20);
end

for ii_n=1:N2D
    
    if isempty(find(Gnode1==ii_n))==1  % exclude the nodes in the punched hole area
        Jt(ii_n)=0;
        Jte(ii_n)=0;
    else
        x_g=Nodex(ii_n);
        y_g=Nodey(ii_n);
        rho_g=sqrt((x_g-L_side/2)^2+(y_g-L_side/2)^2);  %%distance to the center
        %%% the factor of 0.7 below shall be adaptively adjusted
        if rho_g>0.7*L_side/2   % far away from edge, tunneling in z direction 
            Ec=Ec3D(nox+1:nloc+nox+1,ii_n);         %%%band profile in local mesh grid
            Fn=Fn3D(nox+1:Nz,ii_n); 
            Fn_end=Fn(length(Fn));                  %%%quasi fermi level at drain side
            Ec=Ec+phid;
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
            
            %% thermionic emission current over the top of barrier
            N_grid=200;
            EnergyGrid=linspace(0,0.5,N_grid);
            dEnergy=EnergyGrid(2)-EnergyGrid(1);
            TE_bridge=q*Richardson*T/kB*(log((1+exp((-phis-EnergyGrid)*q/kB/T)))-log(1+exp((Fn_end-phis-EnergyGrid)*q/kB/T)));
            Jte(ii_n)=abs(dEnergy*sum(TE_bridge));                           %%thermionic current
            %% end of thermionic emission current over the top of barrier
        else   % finde the tunneling path first, then calculate the tunneling current
            for ii_t=1:length(theta) % iterate over different theta directions for finding the tunneling path
                x_grid=(L_side/2-rho_g)+z*cos(theta(ii_t)); % polar coordingate, polar symmetry assumed
                z_grid=z*sin(theta(ii_t));
                %% Ecz2D(X,Y) is the Ec of a vertical cross section through the hole center
                Ec_nonlocal(:,ii_t)=interp2(X,Y,Ecz2D,x_grid,z_grid);
                Ec_5(ii_t)=Ec_nonlocal(5,ii_t);
            end
            %% the tunneling path determined by minimum Ec at the 5th node
            %% the value of 5 is emprical and adjustable.
            ind_path=find(Ec_5==min(Ec_5));
            Ec=Ec_nonlocal(:,ind_path(length(ind_path)));
            phis=Ec(1);
            for ii_z=1:length(z)
                if ii_z==1
                    Tr(ii_z)=1;
                else
                    %%% Transmission computed by WKB formalism
                    Tr(ii_z)=exp(-1*abs(2/hbar*trapz(z(1:ii_z),abs(sqrt(2*Mass*q*(phis-Ec(1:ii_z)))))));
                end
            end
            dE=Ec(1:length(Ec)-1)-Ec(2:length(Ec));
            T_bridge=q*Richardson*T/kB*(log((1+exp((-Ec)*q/kB/T)))-log(1+exp((Fn_end-Ec)*q/kB/T)));%%%For tunneling
            Jt(ii_n)=abs(sum(Tr(1:length(Tr)-1)'.*dE.*T_bridge(1:length(T_bridge)-1)));       %%tunneling current
            
            %% thermionic emission current over the top of barrier
            N_grid=200;
            EnergyGrid=linspace(0,0.5,N_grid);
            dEnergy=EnergyGrid(2)-EnergyGrid(1);
            TE_bridge=q*Richardson*T/kB*(log((1+exp((-phis-EnergyGrid)*q/kB/T)))-log(1+exp((Fn_end-phis-EnergyGrid)*q/kB/T)));
            Jte(ii_n)=abs(dEnergy*sum(TE_bridge));
            % end of thermionic emission current over the top of barrier
        end
    end
end

%% It seems Wenchao's treatment below omits the thermionic emission current
J_total=sum((Jt)'.*NodeA)/sum(NodeA);
            
            
            