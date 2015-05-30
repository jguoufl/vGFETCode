%%%%%% solve the Graphene based Field Effect Transistor, 3D Poisson-DD
%%%%%% Wenchao Chen, July 2012
%%%%%% DD solved in the thin film channel 
%%%%%% Poisson: solved in the whole channel+oxide region+graphene region
%%%%%% Poisson: charge exists in the channel region, graphene layer
%%% Changed by JG Apr 15, to treat metal instead of graphene source.
%% Ec_bias solved by Poisson.m is the Vacuum energy minus metal (gate, drain) workfunction
%% the zero energy (reference) is defined as the source Fermi level
%% Note by JG (Apr 15): This version works for any tSr for electrostatics (flag_equ==0),
%% but it works only for tSr=0 (one grid point at z=0 as Source) for transport (flag_equ~=0)

clc
clear all
close all

inp;    % set up input parameters
init;   % initialization

Ec3D=zeros(Nz,N2D,Ng_step+1,Nd_step+1);
Ne3D=zeros(Nz,N2D,Ng_step+1,Nd_step+1);
Id=zeros(Ng_step+1,Nd_step+1);

ind_ch=((nox+1)*N2D+1):N3D;  %%index of channel nodes
ind_Sr=[];  % intialization
for ii_Sr=1:Nz_Sr
    ind_Sr=[ind_Sr (nox+ii_Sr-1)*N2D+Gnode1];        %%index of Source nodes
end

% exclude the source nodes from channel
Acommon=intersect(ind_ch,ind_Sr); 
ind_ch=setxor(ind_ch,Acommon); % exclude the source nodes from channel
N_Sr=length(ind_Sr); % number of source contact nodes
Nch=Nz-(nox+1);              %%number of nodes in channel in z direction
phi_gc_plot=sparse(N2D,1);

criterion_outer=1e-3;

for ii_vg=1:(Ng_step+1)
    Vg_bias=Vg0+(ii_vg-1)*Vg_step;
    Evacg=-Vg_bias;
    %%Ecg=phig-Vg_bias;
    Vg_plot(ii_vg)=Vg0+(ii_vg-1)*Vg_step;
    
    for ii_vd=1:(Nd_step+1)
        Vd_bias=Vd0+(ii_vd-1)*Vd_step;
        Vd_plot(ii_vd)=Vd_bias;
        Evacd=-Vd_bias;
        %%Ecd=phid-Vd_bias;
        
        Ne_bias=sparse(N3D,1);
        Fn_bias=sparse(N3D,1);
        
        %%%Set Boundary Conditions For each loop%%%
        %%For Bottom Boundary, Diriclet B.C. For Gate
        S3D(1:N2D,:)=sparse(N2D,N3D);
        K3D(1:N2D,:)=sparse(N2D,N3D);
        K3D(1:N2D,1:N2D)=sparse(diag(ones(1,N2D)));
        %%Bd3D(1:N2D)=Ecg;
        Bd3D(1:N2D)=Evacg;
        %%For Top Boundary, Diriclet B.C. For Drain
        S3D(((Nz-1)*N2D+1):N3D,:)=sparse(N2D,N3D);
        K3D(((Nz-1)*N2D+1):N3D,:)=sparse(N2D,N3D);
        K3D(((Nz-1)*N2D+1):N3D,((Nz-1)*N2D+1):N3D)=sparse(diag(ones(1,N2D)));
        %%Bd3D(((Nz-1)*N2D+1):N3D)=Ecd;
        Bd3D(((Nz-1)*N2D+1):N3D)=Evacd;
        
        %% b.c. for the metal source, JG, Apr 15.
        S3D(ind_Sr,:)=sparse(N_Sr,N3D);
        K3D(ind_Sr,:)=sparse(N_Sr,N3D);
        K3D(ind_Sr,ind_Sr)=sparse(diag(ones(1,N_Sr)));
        %%Bd3D(1:N2D)=Ecg;
        Evacs=0; % source WF minus metal (D&G) WF
        Bd3D(ind_Sr)=Evacs;
        
        %%End of Boundary Conditions%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% the initial guess as Laplace solution
        %[Ec_bias]=poisson(0,0,1,K3D,S3D,Bd3D,N2D,nox,Gnode1,Vs,sz1,sz2,phid);
        Ec_bias=zeros(N3D,1);
        % initial guess for Fn: linear drop over z position, in channel
        % region
        Fn_z=linspace(-Vs, -Vd_bias, Nz-nox-1);

        for ii_fn=1:length(Fn_z)
            N_br=((nox+1+ii_fn-1)*N2D+1):((nox+1+ii_fn+1-1)*N2D);
            Fn_bias(N_br)=Fn_z(ii_fn);
        end
        
        error_outer=1;
        %%%%%%%%%%%%% self-consistent iteration
        while error_outer>criterion_outer
            %%%%% solve the non linear Poisson equation
            Ec_bias_old=Ec_bias;            
            [Ec_bias]=poisson(Fn_bias,Ec_bias_old,0, K3D, S3D,Bd3D,N2D,nox,Gnode1,Vs,sz1,sz2,phid);
            
            if flag_equ==1  % no transport for equilibrium electrostatics
                break
            end
            
            %%% Ec_bias is the vacuum energy minus metal work function (gate & drain)
            
            phi_gc=Ec_bias(ind_Sr)+phid;                  %%shtokky barrier height, graphene-channel, drain workfunction and grpahene dirac point are the same.
            phi_gc_plot(ind_Sr)=Ec_bias(ind_Sr)+phid;
            Fn_bias(ind_Sr+N2D)=0; %% question JG, Apr 15: Why need to "+N2D", shall be removed? 
            %%%%% quasi-Fermi level,solve 3D DD equation
            [Ne_bias(ind_ch),Fn_bias(ind_ch),Jn_f,I_f,Jn_l,I_l]=charge(Ec_bias(ind_ch),Fn_bias(ind_ch),Ele,N2D,NodeA,deltzDD,Nch,Gnode1,phi_gc,phid);
            
            error_outer=max(abs(full(Ec_bias_old-Ec_bias)))
        end
        
        
        Ec3D(:,:,ii_vg,ii_vd)=reshape(full(Ec_bias),N2D,Nz)';  %%write the column first, 
        Ne3D(:,:,ii_vg,ii_vd)=reshape(full(Ne_bias),N2D,Nz)';
        Fn3D(:,:,ii_vg,ii_vd)=reshape(full(Fn_bias),N2D,Nz)';
        Jt(ii_vg,ii_vd)=TCurrent(Ec3D(:,:,ii_vg,ii_vd),Fn3D(:,:,ii_vg,ii_vd),Gnode1,NodeA,nox,N2D,phid,z_final(nox+1:Nz),nloc,Nz);  %%%For D=0
        
        for ii_p=nox+1:Nz
            Ec_p=Ec3D(ii_p,:,ii_vg,ii_vd)+phid;
            vis(:,1)=[Nodex]';      %% in m
            vis(:,2)=[Nodey]';      %% in m
            vis(:,3)=Ec_p;
            [xlin ylin Ec2D]=prz(vis);
            Nx_plot=length(xlin);
            Ny_plot=length(ylin);
            x_plotz(Nx_plot*(ii_p-nox-1)+1:Nx_plot*(ii_p-nox))=xlin;
            z_plotz(Nx_plot*(ii_p-nox-1)+1:Nx_plot*(ii_p-nox))=z_final(ii_p);
            Ec_plotz(Nx_plot*(ii_p-nox-1)+1:Nx_plot*(ii_p-nox))=Ec2D(ceil(Ny_plot/2),:);              %%at certain y
        end
        vis1(:,1)=x_plotz;      %% in m
        vis1(:,2)=z_plotz;      %% in m
        vis1(:,3)=Ec_plotz;  %%carrier density at the beginning of the channel,column vector
        [xlEcz ylEcz Ecz2D]=prz(vis1);    % xlin and ylin in m.
        
%%        [Jn_gt,I_gt(ii_vg,ii_vd)]=chargeGT(Ec_bias(ind_ch),Fn_bias(ind_ch),Ele,N2D,NodeA,deltzDD,Nch,Gnode1,phi_gc,phid,J_tunneling);

        if flag_equ~=1
            [J_total(ii_vg,ii_vd),J_tunneling]=Tpath(Ec3D(:,:,ii_vg,ii_vd),Fn3D(:,:,ii_vg,ii_vd),xlEcz,ylEcz,Ecz2D,Nodex,Nodey,Gnode1,NodeA,nox,N2D,phid,z_final(nox+1:Nz),nloc,Nz,Vg_bias);
            I(ii_vg,ii_vd)=I_l;
            phig_out(ii_vg,ii_vd)=phi_gc(1);
        end
    end % end of the drain voltage loop
end % end of the gate voltage loop
            
%%%%%Visualization%%%%%%%%%
draw;
% I_f
% I_l