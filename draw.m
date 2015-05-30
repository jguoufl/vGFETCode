figure(1)   %% the carrier density of the first sheet in the channel
vis(:,1)=[Nodex]';      %% in m
vis(:,2)=[Nodey]';      %% in m
%vis(:,3)=Fn_bias;
vis(:,3)=Ne3D(nox+2,:,ii_vg,ii_vd);  %%carrier density at the beginning of the channel,column vector
[xlin ylin Ne2D]=pr(vis);    % xlin and ylin in m.
% 
% figure(2)   %% the shttoky barrier height
% vis(:,1)=[Nodex]';      %% in m
% vis(:,2)=[Nodey]';      %% in m
% %vis(:,3)=Fn_bias;
% vis(:,3)=phi_gc_plot;  %%the shttoky barrier height at the beginning of the channel,column vector
% [xlin ylin Phig2D]=pr(vis);    % xlin and ylin in m.
% 
figure(3)   %% Current Density at the first channel sheet
vis(:,1)=[Nodex]';      %% in m
vis(:,2)=[Nodey]';      %% in m
%vis(:,3)=Fn_bias;
vis(:,3)=J_tunneling;  %%Current Density at the first sheet at the beginning of the channel,column vector
[xlin ylin J_tunneling2D]=pr(vis);    % xlin and ylin in m.

% figure(22)   %% Current Density at the first channel sheet
% vis(:,1)=[Nodex]';      %% in m
% vis(:,2)=[Nodey]';      %% in m
% %vis(:,3)=Fn_bias;
% vis(:,3)=Jn_gt;  %%Current Density at the first sheet at the beginning of the channel,column vector
% [xlin ylin Jn_gt2D]=pr(vis);    % xlin and ylin in m.
% 
% figure(4)   %% Current Density at the last channel sheet
% vis(:,1)=[Nodex]';      %% in m
% vis(:,2)=[Nodey]';      %% in m
% %vis(:,3)=Fn_bias;
% vis(:,3)=-Jn_l;  %%Current Density at the last sheet at the beginning of the channel,column vector
% [xlin ylin Jn_l2D]=pr(vis);    % xlin and ylin in m.
figure(5)
semilogy(Vg_plot,J_total,'r','linewidth',[2])
xlabel('Vg (V)','fontsize',[20])
ylabel('Jt (A/m^2)','fontsize',[20])
set(gca,'linewidth',[2],'fontsize',[20])

% figure(11)
% semilogy(Vg_plot,Jt,'r','linewidth',[2])
% xlabel('Vg (V)','fontsize',[20])
% ylabel('J (A/m^2)','fontsize',[20])
% set(gca,'linewidth',[2],'fontsize',[20])

% figure(6)
% plot(Vg_plot,phig_out,'r','linewidth',[2])
% xlabel('Vg (V)','fontsize',[20])
% ylabel('Barrier (eV)','fontsize',[20])
% set(gca,'linewidth',[2],'fontsize',[20])


% Ec_near=Ec3D(nox+1:Nz,149,1,1)+phid;
% Fn_near=Fn3D(nox+1:Nz,149,1,1);
% z_plot=z_final(nox+1:Nz);
% 
% Ec_far=Ec3D(nox+1:Nz,2,1,1)+phid;
% Fn_far=Fn3D(nox+1:Nz,2,1,1);
% 
% figure(8)
% plot(z_plot*1e9,Ec_near,'b','linewidth',[2])
% hold on
% plot(z_plot*1e9,Ec_far,'r','linewidth',[2])
% xlabel('z (nm)','fontsize',[20])
% ylabel('Ec (eV)','fontsize',[20])
% set(gca,'linewidth',[2],'fontsize',[20])
% axis([0,202,-0.5,1]);
oxide_only=0;  % visualize oxide only
if oxide_only==1 
    Nstart=nox+1;    
else
    Nstart=1;    
end
for ii_p=Nstart:Nz
    Ec_p=Ec3D(ii_p,:,Ng_step+1,Nd_step+1)+phid;
    vis(:,1)=[Nodex]';      %% in m
    vis(:,2)=[Nodey]';      %% in m
    vis(:,3)=Ec_p;
    [xlin ylin Ec2D]=prz(vis);
    Nx_plot=length(xlin);
    Ny_plot=length(ylin);
    x_plotz(Nx_plot*(ii_p-Nstart)+1:Nx_plot*(ii_p-Nstart+1))=xlin;
    z_plotz(Nx_plot*(ii_p-Nstart)+1:Nx_plot*(ii_p-Nstart+1))=z_final(ii_p);
    Ec_plotz(Nx_plot*(ii_p-Nstart)+1:Nx_plot*(ii_p-Nstart+1))=Ec2D(ceil(Ny_plot/2),:);              %%at certain y
end
clear vis1
figure(10)   %% the potential profile in x-z plane
vis1(:,1)=x_plotz;      %% in m
vis1(:,2)=z_plotz;      %% in m
%vis(:,3)=Fn_bias;
vis1(:,3)=Ec_plotz;  %%Potential
[xlEcz ylEcz Ecz2D]=pr(vis1);    % xlin and ylin in m.

Ny_cut= sum(ylEcz<0)+1;  % the y index of the graphene layer
figure(12)
Ec_tilt=Ecz2D(Ny_cut,:)';
plot((xlEcz-max(xlEcz)/2)*1e9,Ec_tilt,'linewidth',[2]);
set(gca,'fontsize',[20],'linewidth',[2]);
xlabel('X [nm]');
ylabel('E [eV]');



