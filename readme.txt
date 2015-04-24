init.m: for 3D poisson matrix construction
charge.m: solver for 3D drift diffusion
TCurrent.m: for calculating tunneling current for device without hole
Tpath.m: for calculating tunneling current for device with hole.


input files: 

GFETd0.e and GFETd0.n for device without graphene hole region
GFEThole for device with hole region
GFETholefine for device with hole region, finer grid

%% Note by JG (Apr 15): This version works for any tSr for electrostatics (flag_equ==0),
%% but it works only for tSr=0 (one grid point at z=0 as Source) for transport (flag_equ~=0)