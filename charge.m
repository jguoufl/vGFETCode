% 3D, unipolar drift-diffusion solver
% Wenchao Chen, July 2012
function [Ne_old,Fn_bias,Jn_f,I_f,Jn_l,I_l]=charge(U,Fn_bias_old,Ele,N2D,NodeA,deltzDD,Nch,Gnode1,phi_gc,phid)

%%only channel region is considered, oxide and graphene sheet are excluded.

global epso m0 kBT
global hbar q Ne0 mu

%%charge density for setting boundary conditions%%
Nss=Ne0*fermi(-phi_gc./kBT,1,1/2);
Ndd=Ne0*fermi(-phid/kBT,1,1/2);

%%%%%%%%%%%%%%%End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_n=length(U);      %% the total number of nodes, in channel region
N_e=length(Ele);    %% the total number of elements

%%%%%%% set up the matrices elements for the DD Eqn.
KD=sparse(N_n,N_n);     %%the Laplace operator matrix
BD=sparse(N_n,1);       %%the boundary condition vector

K2D=sparse(N2D,N2D);
S2Dup=sparse(N2D,N2D);  %%up sheet S2D
S2Dr=sparse(N2D,N2D);   %%right S2D
S2Ddown=sparse(N2D,N2D);%%down sheet S2D
S2Dl=sparse(N2D,N2D);   %%left S2D

zetan=(Fn_bias_old-U)./kBT;
deg_fac_c=fermi(zetan,1,1/2)./fermi(zetan,1,-1/2);    %treat degeneratcy and non-parabolic band structure.

for ii_z=2:(Nch-1)            %%the channel region, without the boundary node
    Nos=(ii_z-1)*N2D;         %%node index offset
    zl=(deltzDD(ii_z-1)+deltzDD(ii_z))/2;
    for ii_e=1:N_e
        %%%%%%%% compute the kinetic energy matrix %%%%%%%%%%%
        %%% the 1st edge
        deg_fac=(deg_fac_c(Ele(ii_e).n1+Nos)+deg_fac_c(Ele(ii_e).n2+Nos))/2;
        delt=(U(Ele(ii_e).n2+Nos)-U(Ele(ii_e).n1+Nos))/(kBT*deg_fac);
        c1=mu*(kBT*deg_fac)*Ele(ii_e).d1/Ele(ii_e).l1*Bern(delt);     % coefficient in Schaff-Gummel, for electron flux * (-1)
        c2=mu*(kBT*deg_fac)*Ele(ii_e).d1/Ele(ii_e).l1*Bern(-delt);    % coefficient in Schaff_gunnel 
    
        K2D(Ele(ii_e).n1,Ele(ii_e).n1)=K2D(Ele(ii_e).n1,Ele(ii_e).n1)-c1;
        K2D(Ele(ii_e).n1,Ele(ii_e).n2)=K2D(Ele(ii_e).n1,Ele(ii_e).n2)+c2;
        K2D(Ele(ii_e).n2,Ele(ii_e).n2)=K2D(Ele(ii_e).n2,Ele(ii_e).n2)-c2;
        K2D(Ele(ii_e).n2,Ele(ii_e).n1)=K2D(Ele(ii_e).n2,Ele(ii_e).n1)+c1;
    
    
        %%% the 2nd edge
        deg_fac=(deg_fac_c(Ele(ii_e).n2+Nos)+deg_fac_c(Ele(ii_e).n3+Nos))/2;
        delt=(U(Ele(ii_e).n3+Nos)-U(Ele(ii_e).n2+Nos))/(kBT*deg_fac);    
        c1=mu*(kBT*deg_fac)*Ele(ii_e).d2/Ele(ii_e).l2*Bern(delt);     % coefficient in Schaff-Gummel
        c2=mu*(kBT*deg_fac)*Ele(ii_e).d2/Ele(ii_e).l2*Bern(-delt);    % coefficient in Schaff_gunnel 
        
        K2D(Ele(ii_e).n2,Ele(ii_e).n2)=K2D(Ele(ii_e).n2,Ele(ii_e).n2)-c1;
        K2D(Ele(ii_e).n2,Ele(ii_e).n3)=K2D(Ele(ii_e).n2,Ele(ii_e).n3)+c2;
        K2D(Ele(ii_e).n3,Ele(ii_e).n3)=K2D(Ele(ii_e).n3,Ele(ii_e).n3)-c2;
        K2D(Ele(ii_e).n3,Ele(ii_e).n2)=K2D(Ele(ii_e).n3,Ele(ii_e).n2)+c1;

    
        %%% the 3rd edge
        deg_fac=(deg_fac_c(Ele(ii_e).n3+Nos)+deg_fac_c(Ele(ii_e).n1+Nos))/2;
        delt=(U(Ele(ii_e).n1+Nos)-U(Ele(ii_e).n3+Nos))/(kBT*deg_fac);
        c1=mu*(kBT*deg_fac)*Ele(ii_e).d3/Ele(ii_e).l3*Bern(delt);     % coefficient in Schaff-Gummel
        c2=mu*(kBT*deg_fac)*Ele(ii_e).d3/Ele(ii_e).l3*Bern(-delt);    % coefficient in Schaff_gunnel 
    
        K2D(Ele(ii_e).n3,Ele(ii_e).n3)=K2D(Ele(ii_e).n3,Ele(ii_e).n3)-c1;
        K2D(Ele(ii_e).n3,Ele(ii_e).n1)=K2D(Ele(ii_e).n3,Ele(ii_e).n1)+c2;
        K2D(Ele(ii_e).n1,Ele(ii_e).n1)=K2D(Ele(ii_e).n1,Ele(ii_e).n1)-c2;
        K2D(Ele(ii_e).n1,Ele(ii_e).n3)=K2D(Ele(ii_e).n1,Ele(ii_e).n3)+c1; 
    end
    K2D=zl*K2D;
    for ii_n=1:N2D
        deg_fac=(deg_fac_c(ii_n+Nos)+deg_fac_c(ii_n+Nos+N2D))/2;
        delt=(U(ii_n+Nos+N2D)-U(ii_n+Nos))/(kBT*deg_fac);
        c1=mu*(kBT*deg_fac)*NodeA(ii_n)/deltzDD(ii_z)*Bern(delt);     % coefficient in Schaff-Gummel, for electron flux * (-1)
        c2=mu*(kBT*deg_fac)*NodeA(ii_n)/deltzDD(ii_z)*Bern(-delt);    % coefficient in Schaff_gunnel
        S2Dup(ii_n,ii_n)=-c1;
        S2Dr(ii_n,ii_n)=c2;
        
        deg_fac=(deg_fac_c(ii_n+Nos-N2D)+deg_fac_c(ii_n+Nos))/2;
        delt=(U(ii_n+Nos-N2D)-U(ii_n+Nos))/(kBT*deg_fac);
        c1=mu*(kBT*deg_fac)*NodeA(ii_n)/deltzDD(ii_z-1)*Bern(delt);     % coefficient in Schaff-Gummel, for electron flux * (-1)
        c2=mu*(kBT*deg_fac)*NodeA(ii_n)/deltzDD(ii_z-1)*Bern(-delt);    % coefficient in Schaff_gunnel
        S2Ddown(ii_n,ii_n)=-c1;
        S2Dl(ii_n,ii_n)=c2;
    end
    
    if ii_z==2                                                          %%for the second module of KD            
        S2Dtemp=sparse(N2D,N2D);
        S2Dtemp(Gnode1,Gnode1)=S2Ddown(Gnode1,Gnode1);
        KD(Nos+1:Nos+N2D,Nos+1:Nos+N2D)=K2D+S2Dup+S2Dtemp;
        
        KD(Nos+Gnode1,Nos-N2D+Gnode1)=S2Dl(Gnode1,Gnode1);
        KD(Nos+1:Nos+N2D,Nos+N2D+1:Nos+2*N2D)=S2Dr;
    else
        KD(Nos+1:Nos+N2D,Nos+1:Nos+N2D)=K2D+S2Dup+S2Ddown;
        KD(Nos+1:Nos+N2D,Nos-N2D+1:Nos)=S2Dl;
        KD(Nos+1:Nos+N2D,Nos+N2D+1:Nos+2*N2D)=S2Dr;
    end
    
    S2Dup=sparse(N2D,N2D);                                %%reset the matrix for the next step
    S2Dr=sparse(N2D,N2D);
    S2Ddown=sparse(N2D,N2D); 
    S2Dl=sparse(N2D,N2D); 
    K2D=sparse(N2D,N2D);             
end

% %%%for the first sheet in the channel region; Boundary 
ii_z=1;
Nos=(ii_z-1)*N2D;

    for ii_ng=1:length(Gnode1)                           %%set the terms corresponding to graphene-channel, to 1
        K2D(Gnode1(ii_ng),:)=sparse(1,N2D);
        K2D(Gnode1(ii_ng),Gnode1(ii_ng))=1;
    end
    index=0;
    for ii_n=1:N2D
        index=find(Gnode1==ii_n);
        if isempty(index)
            deg_fac=(deg_fac_c(ii_n+Nos)+deg_fac_c(ii_n+Nos+N2D))/2;
            delt=(U(ii_n+Nos+N2D)-U(ii_n+Nos))/(kBT*deg_fac);
            c1=mu*(kBT*deg_fac)*NodeA(ii_n)/deltzDD(ii_z)*Bern(delt);     % coefficient in Schaff-Gummel, for electron flux * (-1)
            c2=mu*(kBT*deg_fac)*NodeA(ii_n)/deltzDD(ii_z)*Bern(-delt);    % coefficient in Schaff_gunnel
            S2Dup(ii_n,ii_n)=-c1;
            S2Dr(ii_n,ii_n)=c2;
        end
    end
%     S2Dup_temp=sparse(N2D,N2D);
%     S2Dr_temp=sparse(N2D,N2D);
%     S2Dup_temp(Gnode1,Gnode1)=S2Dup(Gnode1,Gnode1);
%     S2Dr_temp(Gnode1,Gnode1)=S2Dr(Gnode1,Gnode1);
    
    KD(Nos+1:Nos+N2D,Nos+1:Nos+N2D)=K2D+S2Dup;
    KD(Nos+1:Nos+N2D,Nos+N2D+1:Nos+2*N2D)=S2Dr;
    BD(Nos+Gnode1)=Nss;         %%%%carrier density at graphene sheet node
%%%% End of bottom boundary%%%%

%%%for the last sheet in the channel region, Boundary
KD(N_n-N2D+1:N_n,N_n-N2D+1:N_n)=sparse(eye(N2D,N2D));
BD(N_n-N2D+1:N_n)=Ndd;      %%%%carrier density at the drain contact
%%% End of top boundary%%%%

%%%Solve for the Electron Density%%%%
Ne_old=KD\BD;

%%%Compute the quasi Fermi levels%%%
Fn_bias=zeros(N_n,1);
zetan=anti_dummy(Ne_old./Ne0,1,1/2);
Fn_bias=kBT*zetan+U;

%%%%%%%%Calculate Id For The First Sheet%%%%%%%%%%
ii_z=1;                                           %%For the first sheet
Nos=(ii_z-1)*N2D;
    for ii_n=1:N2D
        deg_fac=(deg_fac_c(ii_n+Nos)+deg_fac_c(ii_n+Nos+N2D))/2;
        delt=(U(ii_n+Nos+N2D)-U(ii_n+Nos))/(kBT*deg_fac);
        c1_f(ii_n)=mu*(kBT*deg_fac)/deltzDD(ii_z)*Bern(delt);     % coefficient in Schaff-Gummel, for electron flux * (-1)
        c2_f(ii_n)=mu*(kBT*deg_fac)/deltzDD(ii_z)*Bern(-delt);    % coefficient in Schaff_gunnel
%         S2Dup(ii_n,ii_n)=-c1;
%         S2Dr(ii_n,ii_n)=c2;
    end
    
Jn_f=q*(-c1_f'.*Ne_old(Nos+1:Nos+N2D)+c2_f'.*Ne_old(Nos+N2D+1:Nos+2*N2D));      % the electron current density, position dependent, at the first sheet
I_f=full(sum(Jn_f.*NodeA)/sum(NodeA));


%%%%%%%%Calculate Id For The Last Sheet%%%%%%%%%%
ii_z=Nch-1;                                           %%For the Last sheet
Nos=(ii_z-1)*N2D;
    for ii_n=1:N2D
        deg_fac=(deg_fac_c(ii_n+Nos)+deg_fac_c(ii_n+Nos+N2D))/2;
        delt=(U(ii_n+Nos+N2D)-U(ii_n+Nos))/(kBT*deg_fac);
        c1_l(ii_n)=mu*(kBT*deg_fac)/deltzDD(ii_z)*Bern(delt);     % coefficient in Schaff-Gummel, for electron flux * (-1)
        c2_l(ii_n)=mu*(kBT*deg_fac)/deltzDD(ii_z)*Bern(-delt);    % coefficient in Schaff_gunnel
%         S2Dup(ii_n,ii_n)=-c1;
%         S2Dr(ii_n,ii_n)=c2;
    end
    
Jn_l=q*(-c1_l'.*Ne_old(Nos+1:Nos+N2D)+c2_l'.*Ne_old(Nos+N2D+1:Nos+2*N2D));        % the electron current density, position dependent, at the first sheet
I_l=full(sum(Jn_l.*NodeA)/sum(NodeA));

