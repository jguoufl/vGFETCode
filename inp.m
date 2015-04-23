%%%%%% the input parameters for the 3D Poisson-DD solver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of geometric parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Physical constant
global epso m0 kBT
global hbar q Ne0
epso=8.854e-12;
m0=9.11e-31;
kBT=0.0259;
hbar=1.055e-34;
q=1.6e-19;
Ne0=1e22;%2*(m0*kBT*q/(2*pi*hbar^2))^(3/2); %1e27; % the charge density constant for 3D bulk material
mu=1.5e-4; %1.5e-4;            % m^2/Vs

Vs=0;               % the source (Graphene) voltage
Vg0=25;             % the bottom gate voltage
Vg_step=5;
Ng_step=0;

Vd0=2;
Vd_step=0.5;
Nd_step=0;          % the drain (the top electrode) voltage

%%%%%% input the contact work function - the semiconductor affinity
phig=0;           % Gate work function minus semiconductor affinity, flat band bias
phi_dirac=0.88;      % The dirac(graphene)-affinity(channel)
phid=0.88;           % The drain-channel SB height

%Geometric parameters
Lch=200*1e-9;       % organic channel thickness
epso1=6;            % organic channel dielectric contant
tins=20*1e-9;       % bottom insulator thickness, nominal value in m
epsoins=4;          % bottom insulator dielectric constant

%%%% the simulated region and the grid spacing
sz1=0.5*1e-9      % initial grid in channel, in m
szcoa=2*1e-9    % initial grid for non-uniform coarse grid, in channel
sz2=0.5*1e-9      % initial grid in oxide,   in m

alphaz1=1.1;    % the grid ratio in the thin film channel
alphaz2=1.1;    % the grid ratio in the gate insulator

locpath=20e-9;  % local tunneling path length
% end of input parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%