%%%%%%% the initialization of 3D Poisson-DD solver
%%%% generate the non-uniform grid in z direction

%%%%Local Path Region for Tunneling%%%%%
nloc=round(locpath/sz1);
znloc=sz1*(0:nloc-1);

%%Coarse z-gride in channel region
ncoa=round(log(1-(Lch-locpath)*(1-alphaz1)/szcoa)/log(alphaz1));  %%number of grid space in z for coarse grid,node: ncoa+1
for ii=1:(ncoa+1)
    if ii==1
        zncoa(ii)=sz1*nloc;
    else
        zncoa(ii)=sz1*nloc+szcoa*(alphaz1^(ii-1)-1)/(alphaz1-1);  %%coarse grid z coordinates
    end
end
reg_mark5=5*ones(1,length(znloc)+length(zncoa));                  %%region mark in z direction, 5 for channel
%%Coarse z-grid in oxide region
nox=round(log(1-tins*(1-alphaz2)/sz2)/log(alphaz2)); 
for ii=1:nox
    zox(nox+1-ii)=-sz2*(alphaz2^(ii)-1)/(alphaz2-1);              %%oxide grid z coordinates
end
reg_mark6=6*ones(1,length(zox));                                  %%region mark in z direction, 6 for oxide

z_final=[zox znloc zncoa];         
reg_mark=[reg_mark6 reg_mark5];
Nz=length(z_final);                                     %%node number in z direction
for zzi=1:length(z_final)-1
    deltz(zzi)=z_final(zzi+1)-z_final(zzi);             %%separation between nodes in z direction
end
deltzDD=deltz((nox+1+1):(Nz-1));                        %%separation in channel region, for DD calculation, excluding the graphene sheet
%%%%% initialize the matrices for solving 3D Poisson
%gridname=input('Please input grid file name: ','s');
gridname='GFEThole';  %%GFET_d1_fine1
Emtr=load([gridname '.e']);
Nmtr=load([gridname '.n']);

ul=1e-8;                % length unit
N_e=length(Emtr);       % the number of 2D elements
N2D=length(Nmtr);       % the number of nodes in a 2D cross section
N3D=N2D*Nz;             % the total number of nodes

%%%%%%%%% calculate 2D matrices as the first step for calculating 3D matrices
%%% initialize an array to record the position and boundary marker of nodes
Nodex=Nmtr(:,2)*ul;  % the x position of the nodes
Nodey=Nmtr(:,3)*ul;  % the y position of the nodes

%%% initialize an array to record the node number, position, area at the
%%% corners of each element
NodeA=zeros(N2D,1);

ii_Marker=1;
for ii_e=1:N_e
    Ele(ii_e).n1=Emtr(ii_e,2)+1;  % the no. of the 1st node starting from 1
    Ele(ii_e).n2=Emtr(ii_e,3)+1;  % the no. of the 2nd node starting from 1
    Ele(ii_e).n3=Emtr(ii_e,4)+1;  % the no. of the 3rd node starting from 1
    
    Ele(ii_e).x1=Nodex(Ele(ii_e).n1);    % the x position of the first node
    Ele(ii_e).y1=Nodey(Ele(ii_e).n1);    % the y position of the first node
    Ele(ii_e).x2=Nodex(Ele(ii_e).n2);    % the x position of the first node
    Ele(ii_e).y2=Nodey(Ele(ii_e).n2);    % the y position of the first node
    Ele(ii_e).x3=Nodex(Ele(ii_e).n3);    % the x position of the first node
    Ele(ii_e).y3=Nodey(Ele(ii_e).n3);    % the y position of the first node
    Ele(ii_e).A=0.5*det([1 Ele(ii_e).x1 Ele(ii_e).y1; 1 Ele(ii_e).x2 Ele(ii_e).y2; 1 Ele(ii_e).x3 Ele(ii_e).y3]); % the area

    if Emtr(ii_e,13)==1        % graphene area, 
        Gnode(3*ii_Marker-2:3*ii_Marker)=[Emtr(ii_e,2)+1 Emtr(ii_e,3)+1 Emtr(ii_e,4)+1];  %%%Mark the node in graphene area
        ii_Marker=ii_Marker+1; %% increase marker for the next step
    end
    %%%% the length of the 3 edges of the triangular element
    Ele(ii_e).l1=sqrt((Ele(ii_e).x1-Ele(ii_e).x2)^2+(Ele(ii_e).y1-Ele(ii_e).y2)^2); 
    Ele(ii_e).l2=sqrt((Ele(ii_e).x2-Ele(ii_e).x3)^2+(Ele(ii_e).y2-Ele(ii_e).y3)^2);
    Ele(ii_e).l3=sqrt((Ele(ii_e).x3-Ele(ii_e).x1)^2+(Ele(ii_e).y3-Ele(ii_e).y1)^2);
    %%%% compute the distances from the circumcenter to the edges
    xctmp=ul*Emtr(ii_e,11);    yctmp=ul*Emtr(ii_e,12);  % the circumcenter position
    Ele(ii_e).d1=sqrt((xctmp-(Ele(ii_e).x1+Ele(ii_e).x2)/2)^2+(yctmp-(Ele(ii_e).y1+Ele(ii_e).y2)/2)^2); 
    Ele(ii_e).d2=sqrt((xctmp-(Ele(ii_e).x2+Ele(ii_e).x3)/2)^2+(yctmp-(Ele(ii_e).y2+Ele(ii_e).y3)/2)^2); 
    Ele(ii_e).d3=sqrt((xctmp-(Ele(ii_e).x3+Ele(ii_e).x1)/2)^2+(yctmp-(Ele(ii_e).y3+Ele(ii_e).y1)/2)^2); 
    %%%% Compute the area around each node
    NodeA(Ele(ii_e).n1)=NodeA(Ele(ii_e).n1)+(Ele(ii_e).l1*Ele(ii_e).d1+Ele(ii_e).l3*Ele(ii_e).d3)/4;
    NodeA(Ele(ii_e).n2)=NodeA(Ele(ii_e).n2)+(Ele(ii_e).l1*Ele(ii_e).d1+Ele(ii_e).l2*Ele(ii_e).d2)/4;
    NodeA(Ele(ii_e).n3)=NodeA(Ele(ii_e).n3)+(Ele(ii_e).l2*Ele(ii_e).d2+Ele(ii_e).l3*Ele(ii_e).d3)/4;
end

Gnode1=unique(Gnode);   %%eleminate node with same index,node index in graphene-channel interface

S2D=sparse(N2D,N2D);    % the 2D overlap matrix 
Kp2D=sparse(N2D,N2D);   % the Laplace operator matrix without dielectric constant for the cross section

for ii_e=1:N_e     
    %%% compute the 2D overlap matrix (x, y) 
    S2D(Ele(ii_e).n1,Ele(ii_e).n1)=S2D(Ele(ii_e).n1,Ele(ii_e).n1)+Ele(ii_e).A/6; %%
    S2D(Ele(ii_e).n1,Ele(ii_e).n2)=S2D(Ele(ii_e).n1,Ele(ii_e).n2)+Ele(ii_e).A/12;
    S2D(Ele(ii_e).n1,Ele(ii_e).n3)=S2D(Ele(ii_e).n1,Ele(ii_e).n3)+Ele(ii_e).A/12;
    S2D(Ele(ii_e).n2,Ele(ii_e).n1)=S2D(Ele(ii_e).n2,Ele(ii_e).n1)+Ele(ii_e).A/12;
    S2D(Ele(ii_e).n2,Ele(ii_e).n2)=S2D(Ele(ii_e).n2,Ele(ii_e).n2)+Ele(ii_e).A/6;
    S2D(Ele(ii_e).n2,Ele(ii_e).n3)=S2D(Ele(ii_e).n2,Ele(ii_e).n3)+Ele(ii_e).A/12;
    S2D(Ele(ii_e).n3,Ele(ii_e).n1)=S2D(Ele(ii_e).n3,Ele(ii_e).n1)+Ele(ii_e).A/12; 
    S2D(Ele(ii_e).n3,Ele(ii_e).n2)=S2D(Ele(ii_e).n3,Ele(ii_e).n2)+Ele(ii_e).A/12; 
    S2D(Ele(ii_e).n3,Ele(ii_e).n3)=S2D(Ele(ii_e).n3,Ele(ii_e).n3)+Ele(ii_e).A/6; 

    
    %% compute the 2D (x, y) Laplace matrix for Poisson %% Kp=-grad(er*grad)
    Kp2D(Ele(ii_e).n1,Ele(ii_e).n1)=Kp2D(Ele(ii_e).n1,Ele(ii_e).n1)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x2-Ele(ii_e).x3)^2+...
            (Ele(ii_e).y2-Ele(ii_e).y3)^2);
    Kp2D(Ele(ii_e).n1,Ele(ii_e).n2)=Kp2D(Ele(ii_e).n1,Ele(ii_e).n2)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x3-Ele(ii_e).x2)*...
            (Ele(ii_e).x1-Ele(ii_e).x3)+(Ele(ii_e).y2-Ele(ii_e).y3)*(Ele(ii_e).y3-Ele(ii_e).y1));
    Kp2D(Ele(ii_e).n1,Ele(ii_e).n3)=Kp2D(Ele(ii_e).n1,Ele(ii_e).n3)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x3-Ele(ii_e).x2)*...
            (Ele(ii_e).x2-Ele(ii_e).x1)+(Ele(ii_e).y2-Ele(ii_e).y3)*(Ele(ii_e).y1-Ele(ii_e).y2));
    Kp2D(Ele(ii_e).n2,Ele(ii_e).n1)=Kp2D(Ele(ii_e).n2,Ele(ii_e).n1)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x3-Ele(ii_e).x2)*...
            (Ele(ii_e).x1-Ele(ii_e).x3)+(Ele(ii_e).y2-Ele(ii_e).y3)*(Ele(ii_e).y3-Ele(ii_e).y1));
    Kp2D(Ele(ii_e).n2,Ele(ii_e).n2)=Kp2D(Ele(ii_e).n2,Ele(ii_e).n2)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x1-Ele(ii_e).x3)^2+...
            (Ele(ii_e).y1-Ele(ii_e).y3)^2);
    Kp2D(Ele(ii_e).n2,Ele(ii_e).n3)=Kp2D(Ele(ii_e).n2,Ele(ii_e).n3)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x1-Ele(ii_e).x3)*...
            (Ele(ii_e).x2-Ele(ii_e).x1)+(Ele(ii_e).y3-Ele(ii_e).y1)*(Ele(ii_e).y1-Ele(ii_e).y2));
    Kp2D(Ele(ii_e).n3,Ele(ii_e).n1)=Kp2D(Ele(ii_e).n3,Ele(ii_e).n1)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x3-Ele(ii_e).x2)*...
            (Ele(ii_e).x2-Ele(ii_e).x1)+(Ele(ii_e).y2-Ele(ii_e).y3)*(Ele(ii_e).y1-Ele(ii_e).y2));
    Kp2D(Ele(ii_e).n3,Ele(ii_e).n2)=Kp2D(Ele(ii_e).n3,Ele(ii_e).n2)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x1-Ele(ii_e).x3)*...
            (Ele(ii_e).x2-Ele(ii_e).x1)+(Ele(ii_e).y3-Ele(ii_e).y1)*(Ele(ii_e).y1-Ele(ii_e).y2));
    Kp2D(Ele(ii_e).n3,Ele(ii_e).n3)=Kp2D(Ele(ii_e).n3,Ele(ii_e).n3)+(1/4/Ele(ii_e).A)*((Ele(ii_e).x1-Ele(ii_e).x2)^2+...
            (Ele(ii_e).y1-Ele(ii_e).y2)^2);   
    
end

%%%the 3D matrices for the 3D Poisson solver %%%%
S3D=sparse(N3D,N3D);
K3D=sparse(N3D,N3D);
Bd3D=sparse(N3D,1);
%% K: (gradient in x,y)*(overlap in z)+(overlap in x,y)*gradient in (z)
Kdiag=Kp2D;  % Laplace operator for Poisson
Sdiag=S2D;   % overlap matrix

%% obtain vectors for the indices and values of non zero elements
[Kdi,Kdj,Kdv]=find(Kdiag);  
[Sdi,Sdj,Sdv]=find(Sdiag);
ldiag=length(Kdi); 

ind_s=1;
%%%%%% construct the 3D S and K in a efficient way using sparse matrix
%%%%%% construction
for ii_z=1:Nz
    if ii_z>nox
        er=epso1;
    else
        er=epsoins;
    end
                       
    if ii_z==1  
        a=deltz(ii_z);
        %%% the upper diagonal block
        S3Di(ind_s:(ind_s+ldiag-1))=Sdi+(ii_z-1)*N2D;
        S3Dj(ind_s:(ind_s+ldiag-1))=Sdj+(ii_z-1)*N2D;
        S3Dv(ind_s:(ind_s+ldiag-1))=a/3*Sdv;        
        K3Di(ind_s:(ind_s+ldiag-1))=Kdi+(ii_z-1)*N2D;
        K3Dj(ind_s:(ind_s+ldiag-1))=Kdj+(ii_z-1)*N2D;
        K3Dv(ind_s:(ind_s+ldiag-1))=er*(a/3*Kdv+1/a*Sdv);
        ind_s=ind_s+ldiag;

        %%% The upper off diagonal block
        S3Di(ind_s:(ind_s+ldiag-1))=Sdi+(ii_z-1)*N2D;
        S3Dj(ind_s:(ind_s+ldiag-1))=Sdj+ii_z*N2D;
        S3Dv(ind_s:(ind_s+ldiag-1))=a/6*Sdv;
        K3Di(ind_s:(ind_s+ldiag-1))=Kdi+(ii_z-1)*N2D;
        K3Dj(ind_s:(ind_s+ldiag-1))=Kdj+ii_z*N2D;
        K3Dv(ind_s:(ind_s+ldiag-1))=er*(a/6*Kdv-1/a*Sdv);
        ind_s=ind_s+ldiag;
        %%% The lower off diagonal block
        S3Di(ind_s:(ind_s+ldiag-1))=Sdi+ii_z*N2D;
        S3Dj(ind_s:(ind_s+ldiag-1))=Sdj+(ii_z-1)*N2D;
        S3Dv(ind_s:(ind_s+ldiag-1))=a/6*Sdv;
        K3Di(ind_s:(ind_s+ldiag-1))=Kdi+ii_z*N2D;
        K3Dj(ind_s:(ind_s+ldiag-1))=Kdj+(ii_z-1)*N2D;
        K3Dv(ind_s:(ind_s+ldiag-1))=er*(a/6*Kdv-1/a*Sdv);      
        ind_s=ind_s+ldiag;
        %%% increase the index
        
    elseif ii_z==Nz
        a=deltz(ii_z-1);
        %%% the diagonal block
        S3Di(ind_s:(ind_s+ldiag-1))=Sdi+(ii_z-1)*N2D;
        S3Dj(ind_s:(ind_s+ldiag-1))=Sdj+(ii_z-1)*N2D;
        S3Dv(ind_s:(ind_s+ldiag-1))=a/3*Sdv;
        K3Di(ind_s:(ind_s+ldiag-1))=Kdi+(ii_z-1)*N2D;
        K3Dj(ind_s:(ind_s+ldiag-1))=Kdj+(ii_z-1)*N2D;
        K3Dv(ind_s:(ind_s+ldiag-1))=er*(a/3*Kdv+1/a*Sdv);
        ind_s=ind_s+ldiag;
   
    else
        a=deltz(ii_z-1);                                   %%seperation in the ii_z grid
        b=deltz(ii_z);                                     %%seperation in the next grid
        %%% the diagonal block
        S3Di(ind_s:(ind_s+ldiag-1))=Sdi+(ii_z-1)*N2D;
        S3Dj(ind_s:(ind_s+ldiag-1))=Sdj+(ii_z-1)*N2D;
        S3Dv(ind_s:(ind_s+ldiag-1))=(a/3+b/3)*Sdv;        
        K3Di(ind_s:(ind_s+ldiag-1))=Kdi+(ii_z-1)*N2D;
        K3Dj(ind_s:(ind_s+ldiag-1))=Kdj+(ii_z-1)*N2D;
        K3Dv(ind_s:(ind_s+ldiag-1))=er*((a/3+b/3)*Kdv+(1/a+1/b)*Sdv);
        ind_s=ind_s+ldiag;

        %%% The upper off diagonal block
        S3Di(ind_s:(ind_s+ldiag-1))=Sdi+(ii_z-1)*N2D;
        S3Dj(ind_s:(ind_s+ldiag-1))=Sdj+ii_z*N2D;
        S3Dv(ind_s:(ind_s+ldiag-1))=b/6*Sdv;
        K3Di(ind_s:(ind_s+ldiag-1))=Kdi+(ii_z-1)*N2D;
        K3Dj(ind_s:(ind_s+ldiag-1))=Kdj+ii_z*N2D;
        K3Dv(ind_s:(ind_s+ldiag-1))=er*(b/6*Kdv-1/b*Sdv);
        ind_s=ind_s+ldiag;
        
        %%% The lower off diagonal block
        S3Di(ind_s:(ind_s+ldiag-1))=Sdi+ii_z*N2D;
        S3Dj(ind_s:(ind_s+ldiag-1))=Sdj+(ii_z-1)*N2D;
        S3Dv(ind_s:(ind_s+ldiag-1))=b/6*Sdv;
        K3Di(ind_s:(ind_s+ldiag-1))=Kdi+ii_z*N2D;
        K3Dj(ind_s:(ind_s+ldiag-1))=Kdj+(ii_z-1)*N2D;
        K3Dv(ind_s:(ind_s+ldiag-1))=er*(b/6*Kdv-1/b*Sdv);  
        ind_s=ind_s+ldiag; 
        %%% increase the index
              
    end
end
S3D=sparse(S3Di,S3Dj,S3Dv);
K3D=sparse(K3Di,K3Dj,K3Dv);
K3D=-1*K3D;                      %%For FEM method, minus sign for dynamic matrix

%% set the boundary conditions

% %%For Bottom Boundary, Diriclet B.C. For Gate
% S3D(1:N2D,:)=sparse(N2D,N3D);
% 
% K3D(1:N2D,:)=sparse(N2D,N3D);
% K3D(1:N2D,1:N2D)=sparse(diag(ones(1,N2D)));
% Bd3D(1:N2D)=-Vg0;
% 
% %%For Top Boundary, Diriclet B.C. For Drain
% S3D(((Nz-1)*N2D+1):N3D,:)=sparse(N2D,N3D);
% 
% K3D(((Nz-1)*N2D+1):N3D,:)=sparse(N2D,N3D);
% K3D(((Nz-1)*N2D+1):N3D,((Nz-1)*N2D+1):N3D)=sparse(diag(ones(1,N2D)));
% Bd3D(((Nz-1)*N2D+1):N3D)=-Vd0+phid;