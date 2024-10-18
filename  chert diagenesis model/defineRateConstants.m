function [rateConstants] = defineRateConstants(A,Ea,R,temp_vec,a,detMin)

% PJF 2021: Define the rate constants for dissolution of each silica phase
% at each depth in the sediment column. Includes the increase in quartz
% dissolution/precipitation as a function of detrital mineral content.

% These are dissolution rate constants; precipitation rate constants used
% later in the model are defined from the from principle of microscopic
% reversibility (i.e. Ksp = k_diss/k_precip, see e.g. Palandri and Kharaka
% 2004.)

    k1 =  A(1) .* exp(-(Ea(1))./(R.*(temp_vec)));                             
    k2 =  A(2) .* exp(-(Ea(2))./(R.*(temp_vec))); 
    k3 =  A(3) .* exp(-(Ea(3))./(R.*(temp_vec))); 
    
% Use an expression of the form log10(k/k0) = ax^b. At 100% det min, rate
% is at maximum, and at 0% det min, rate is about 3.8 orders of magnitude
% less (this comes from the calibration exercise perfomed earlier, see
% figure S10 in manuscript).

suppressionByDetritalMinerals = 10^(a*detMin) / 10^a;  

% This is necessary to make the opal-A persist long enough in the sediment
% column; experimental rate constants are slightly too fast. Makes sense
% in the context of early silica diaenesis work (e.g. Dixit et al, van
% Beusekom et al, Van Cappellen and Qui, etc).
suppressionfromYanchilinaCalibration = 0.05;       

k3 = k3 * suppressionByDetritalMinerals;
k1 = k1 * suppressionfromYanchilinaCalibration;

    
    rateConstants.k1 = k1; %in 1/kyr
    rateConstants.k2 = k2;
    rateConstants.k3 = k3; 
