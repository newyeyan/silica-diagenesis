function [z_real,z,thickness,sed_amt] = defineDepthDomain(dz,N,c,por_0,dens)

% PJF 04-01-2018. Script to define depth intervals for a series of sediment
% boxes with progressively more compacted sediment. Porosity decreases
% according to an exponential function (e.g. Athy 1930, referenced in
% LaRowe et al. 2017). So that there is the same amount of sediment in each
% box, a solver routine is used (there is no exact analytical solution to
% the integral of this exponential function). The depth intervals are then
% used to define the integral, while using the same quantity of sediment in
% a box means the advection becomes simpler.

% por_z = por_0 * exp(-cz), where por = porosity, c = the sediment compaction length scale (m^-1), and z = depth (m)

dz = dz * 1000;             % convert from km to m
z = (0 : dz : N*dz) + dz/2; % the uncorrected depth domain

Vol_ref = dz - por_0/c + (por_0*exp(-c*dz))/c; % Reference volume of sediment in surface (uppermost) box (this is the definite integral between z = 0 and z = 5

opts = optimoptions('fsolve');
opts.Display = 'off';

z_real = zeros(1,length(z)+1);
thickness = zeros(1,length(z));
sed_vol = zeros(1,length(z));


z_real(1) = 0;

for i = 1:length(z)
    fun = @(z_r) abs(z_r - por_0/c + (por_0*exp(-c*z_r))/c - Vol_ref*i);
    z_real(i+1) = fsolve(fun,dz*i,opts);
    thickness(i) = z_real(i+1) - z_real(i);
    sed_vol(i) = z_real(i+1) + (por_0*exp(-c*z_real(i+1)))/c - z_real(i) - (por_0*exp(-c*z_real(i)))/c;
end

sed_amt = sed_vol * dens;


