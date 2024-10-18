function [massStruct,ratioStruct] = preallocateVariables(z,R_SMOW,d18_w_i,overlying_water_Si_conc,porespace,density_water,SolQtz,R_NBS_2928,R_NBS_3028,d30_w_i,beta)
 
% revise from PJF 2021.

% ################### Initialise and preallocate the variables in the mass-balance model  ###################
    % General formulation: There are three solid O bearing phases (opalA, opalCT, quartz), and
    % two dissolved phases (water + Si).  These are both treated separately for 16O
    % and 18O, 28Si and 30Si)

    % Create two structures for ease of use that hold all the variables:
    % 'massStruct' holds the masses of 16O, 18O, 28Si, 30Si in each phase.
    % 'ratioStruct' holds ratios and delta values in each phase.
    
    % First for 16O (solids) :
    opalA_16 = ones(size(z)) * 1e4 ;  % 1e3  mass vectors, units moles/layers
    opalAp1_16 = ones(size(z)) * 0.0;
    
    opalCT_16 = ones(size(z)) *  1e-6  ;             
    opalCTp1_16 = ones(size(z)) * 0.0;             % ranges from d(101) = 4.11 (beginning) to d(101) = 4.055.
    
    qtz_16 = ones(size(z)) *  1e-6; % 1e-3 1   1e3   1e4           % NB: need to include small amount of opalCT, quartz initially to avoid numerical issues.
    qtzp1_16 = ones(size(z)) * 0.0;

   % Then for 18O (solids) :
    opalA_18 = opalA_16 * R_SMOW;  % Give everything the same initial isotope composition (i.e. 0 permil), Mass balance effects are negligble. 
    opalAp1_18 = opalAp1_16 * R_SMOW;
    
    opalCT_18 = opalCT_16 * R_SMOW;
    opalCTp1_18 = opalCTp1_16 * R_SMOW;
    
    qtz_18 = qtz_16 * R_SMOW;
    qtzp1_18 = qtzp1_16 * R_SMOW;

    % Then for the pore fluids
    R0_porewater_O = (d18_w_i/1000 + 1) * R_SMOW;    % Initially all pore water has same isotope compostion, equivalent to a typical (modern) ocean value of -1.1 permil (expressed here as an isotope ratio)  
    tot_water = (porespace * density_water / 18.02);  % How many moles of water per box, a crass approximation. 18.02 = Mr water, density of water in units of g/m3, porespace in units of m3.
    
    porewater_16 = tot_water./(R0_porewater_O+1);      % moles of the individual isotopes per box.
    porewater_18 = tot_water - porewater_16;
    
     porewater_Si = porespace .* 1000 .* SolQtz;                 %  i.e. moles per box.
     porewater_Si_conc = porewater_Si./(porespace*1000);         % units of mol/litre
    porewater_Si_conc(1)= overlying_water_Si_conc;
    porewater_Si(1)= porewater_Si_conc(1) *1000;

    R0_porewater_Si = (d30_w_i/1000 + 1) * R_NBS_3028;  % 30/28 ratio of porewater, also equivalent to a typical (modern) ocean value of 1 permil
    R1_porewater_Si = (d30_w_i/1000 + 1)^beta * R_NBS_2928;  % 29/28 ratio of porewater

    porewater_28 = porewater_Si./(R0_porewater_Si+ 1 +R1_porewater_Si);          % moles of the individual isotopes per box.
    porewater_30 = porewater_28 * R0_porewater_Si ;
    
    massStruct = struct;

    massStruct.porewater_Si = porewater_Si;
    
    massStruct.porewater_16O = porewater_16;
    massStruct.porewater_18O = porewater_18;
    massStruct.opalA_16O = opalA_16;
    massStruct.opalA_18O = opalA_18;
    massStruct.opalCT_16O = opalCT_16;
    massStruct.opalCT_18O = opalCT_18;
    massStruct.quartz_16O = qtz_16;
    massStruct.quartz_18O = qtz_18;

    massStruct.opalA_Si = (massStruct.opalA_16O + massStruct.opalA_18O) / 2;        % Si masses are just defined from the oxygen masses, divided by two (for SiO2 stochiometry). Neglects the influence of structural water.
    massStruct.opalCT_Si = (massStruct.opalCT_16O + massStruct.opalCT_18O) / 2;
    massStruct.quartz_Si = (massStruct.quartz_18O + massStruct.quartz_16O) / 2;

    massStruct.porewater_28Si = porewater_28;
    massStruct.porewater_30Si = porewater_30;
    
    massStruct.opalA_28Si = (massStruct.opalA_16O + massStruct.opalA_18O) / 2 * 1 /(1+R_NBS_3028+R_NBS_2928);
    massStruct.opalA_30Si = (massStruct.opalA_16O + massStruct.opalA_18O) / 2 * R_NBS_3028 /(1+R_NBS_3028+R_NBS_2928);
    massStruct.opalCT_28Si = (massStruct.opalCT_16O + massStruct.opalCT_18O) / 2 * 1 /(1+R_NBS_3028+R_NBS_2928);
    massStruct.opalCT_30Si = (massStruct.opalCT_16O + massStruct.opalCT_18O) / 2 * R_NBS_3028 /(1+R_NBS_3028+R_NBS_2928);
    massStruct.quartz_28Si = (massStruct.quartz_18O + massStruct.quartz_16O) / 2 * 1 /(1+R_NBS_3028+R_NBS_2928);
    massStruct.quartz_30Si = (massStruct.quartz_18O + massStruct.quartz_16O) / 2 * R_NBS_3028 /(1+R_NBS_3028+R_NBS_2928);
  
    
    massStruct.tot_water = tot_water;
    
    ratioStruct = struct;

    ratioStruct.porewater_O = massStruct.porewater_18O./massStruct.porewater_16O;
    ratioStruct.opalA_O = massStruct.opalA_18O./massStruct.opalA_16O;
    ratioStruct.opalCT_O = massStruct.opalCT_18O./massStruct.opalCT_16O;
    ratioStruct.quartz_O = massStruct.quartz_18O./massStruct.quartz_16O;
    ratioStruct.porewater_d18O = (ratioStruct.porewater_O ./ R_SMOW - 1) * 1000;
    ratioStruct.opalA_d18O = (ratioStruct.opalA_O ./ R_SMOW - 1) * 1000;
    ratioStruct.opalCT_d18O = (ratioStruct.opalCT_O ./ R_SMOW - 1) * 1000;
    ratioStruct.quartz_d18O = (ratioStruct.quartz_O ./ R_SMOW - 1) * 1000;

    ratioStruct.porewater_Si = massStruct.porewater_30Si./massStruct.porewater_28Si;
    ratioStruct.opalA_Si = massStruct.opalA_30Si./massStruct.opalA_28Si;
    ratioStruct.opalCT_Si = massStruct.opalCT_30Si./massStruct.opalCT_28Si;
    ratioStruct.quartz_Si = massStruct.quartz_30Si./massStruct.quartz_28Si;
    ratioStruct.porewater_d30Si = (ratioStruct.porewater_Si ./ R_NBS_3028 - 1) * 1000;
    ratioStruct.opalA_d30Si = (ratioStruct.opalA_Si ./ R_NBS_3028 - 1) * 1000;
    ratioStruct.opalCT_d30Si = (ratioStruct.opalCT_Si ./ R_NBS_3028 - 1) * 1000;
    ratioStruct.quartz_d30Si = (ratioStruct.quartz_Si ./ R_NBS_3028 - 1) * 1000;



    
%     preallocatedVariables = 'preallocatedVariables.mat';
%     save(preallocatedVariables)