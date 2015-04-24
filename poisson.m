%%%%%%%%% A general purpose 3D non-linear Poisson equation
%%%%%%%%%%%%%%%%%%%%%%by Wenchao Chen%%%%%%%%%%%%%%%%%%%%%
function [Ec_bias]=poisson(Fn,Ec_old,Laplace_flag, K, S,Bc,N2D,nox,Gnode1,Vs,sz1,sz2,phid)

%%%%%%% Input: 
%% Fn: quasi-Fermi level
%% Ec_old: the old band profile
%% Laplace_flag: 1 for Laplace solution, otherwise for Poisson solution
%% K: the Laplace operator on the R.H.S. of the equation 
%% S: The overlap matrix
%% Bc: the boundary node submatrix on the L.H.S. of the equation
%% N2D: Number of nodes for each layer, nox+1:channel index
%% nox: the last oxide node index in z direction,
%% Gnode1: Graphene node index array. Vs:sorce voltage, determine Ef in
%% Graphene; sz1,sz2:gride in z direction adjacent to graphene layer
%%%%%%% Output:
%% Ec_bias: the band profile
%% Ne_bias: the charge density at equil. computed by
%% Ne0*fermi((Fn-Ec)/kBT,1,1/2),
global epso m0 kBT
global hbar q Ne0
global ind_ch
vF=9.3e5;                      %%Fermi velocity of graphene  
nu_tot=length(K);
 
% Ne_g=sparse(N2D,1);            %%charge density for graphene layer nodes
% Nz_g=nox+1;                     %%the graphene sheet node, positioned at z=0;
% Nz2=Nz_g+1;                     %%start of channel node
% ind_ch=((Nz2-1)*N2D+1):nu_tot;  %%index of channel nodes

%%%%%%%%%%%%%%%For graphene sheet charge density%%%%%%%%%%%
% ind_Sr=(Nz_g-1)*N2D+Gnode1;      %%index of source nodes
%%%%%%%%%%%%%%%%%%%%%End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Laplace_flag==1
    %%%%%% initial guess as Laplace solution
    Ec_bias=real(K\Bc);     % a direct method to solve KX=Bc
else
    error_inner=1; 
    criterion_inner=1e-3;
    Ec_q=phid+Ec_old(ind_ch);         %%with respect to vacum level.
    while error_inner>criterion_inner
        dummy_ro=zeros(nu_tot,1); dummy_ro_prime=zeros(nu_tot,1);
        
        %%% compute the total charge density (electron as positive)
        dummy_ro(ind_ch)=Ne0*fermi((Fn(ind_ch)-Ec_q)./kBT,1,1/2);                 %%for Channel region
        dummy_ro_prime(ind_ch)=-(Ne0/kBT)*fermi((Fn(ind_ch)-Ec_q)./kBT,1,-1/2);                    %%negative value
        
        %%% treatment of charge in graphene source, commented out by JG for
        %%% treating metal contacts
%         Ne_g(Gnode1)=q^2/(pi*hbar^2*vF^2)*sign(-Ec_old(ind_g)-Vs).*(-Vs-Ec_old(ind_g)).^2/(sz1/2+sz2/2);  %%graphene sheet,positive for eletron, negative for hole
%         dummy_ro(ind_g)=Ne_g(Gnode1);
        %%% compute the derivitive of the toatl charge density.
%         dummy_ro_prime(ind_g)=-abs(2*q^2/(pi*hbar^2*vF^2).*abs(-Vs-Ec_old(ind_g))/(sz1/2+sz2/2));  %%negative value, always

        %%%%% Newton-Ralphson solution
        Res=K*Ec_old+(q/epso)*S*(dummy_ro)-Bc;
        Jm=K+(q/epso)*S*spdiags(dummy_ro_prime,0,nu_tot,nu_tot);
        delt_U=-sparse(Jm)\Res;

        %% bound delt_U to eliminate non-physical values
        for i_node=1:length(delt_U)
            if(abs(delt_U(i_node))<=1)
                delt_U(i_node)=delt_U(i_node);
            elseif(1<abs(delt_U(i_node)) & abs(delt_U(i_node)) <3.7)
                delt_U(i_node)=sign(delt_U(i_node))*power(abs(delt_U(i_node)),1/5);
            elseif(abs(delt_U(i_node))>=3.7)
                delt_U(i_node)=sign(delt_U(i_node))*log(abs(delt_U(i_node)));
            end
        end
        
        %%% update the parameters
        error_inner=max(abs(full(delt_U)))
        Ec_bias=Ec_old+delt_U;
        Ec_old=Ec_bias;        
    end
%     Ne_bias=dummy_ro;  %%is necessary? or not?
end

