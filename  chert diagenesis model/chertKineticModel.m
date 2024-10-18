%% silicon and oxygen isotopes during diagenesis, revise from PJF 2022
clear, clc, close all

%% Define the parameter space to be tested

depSetting = 'margin';                  % Depositional setting, can be one of 'margin', 'shelf', or 'abyss' (see LaRowe et al. 2017 for associated parameterisation of sediment thermal conductivty, porosity, compaction, etc.)
overlying_water_Si_conc = 1500 * 10^-6;  % bottom water Si concentration in mol/litre
v_surf = 0.01;                          % surface sedimentation rate, m/kyr .
Sed_Si_pc = 20;                         % percent sediment drymass that is SiO2. opal A
bwT =  25 + 273.15 ;                    % Temperature of bottom waters overlying sediment, C into Kelvin
heatFlow = 0.2;                        % in W/m2, typical values today of the order 0.1 see e.g. Sclater et al. 1980)
mineralContent = 20; % detrital mineral contents to be tested, as percent.
fractionation_si = 1.0021; % fractionation of initial opal A sediment

%% Define fixed and variable parameters. 
% ################### Define basic physical variables ###################

cd(fileparts(matlab.desktop.editor.getActiveFilename)); % Change to folder containing this script (if not there already).

switch depSetting
    case 'shelf'
        c = 0.5 * 10^-3;                                % Sediment compaction length scale, from LaRowe et al. 2017
        por_0 = 0.45;                                   % porosity at sediment-water interface (NB, not the very upper cm). Indicative value from Hantschel and Kauerauf 2009, cited in LaRowe et al. 2017
        Cond_sed = 3.2;                                 % W m^-1 K^-1, typical value from LaRowe et al 2017
        waterColumnDepth = 100;                         % Abritrary water column depth (has very minor impact on solubilities via pressure)
    case 'margin'
        c = 1.7 * 10^-3;                                % Sediment compaction length scale, from LaRowe et al. 2017
        por_0 = 0.74;                                   % porosity at sediment-water interface (NB, not the very upper cm). Indicative value from Hantschel and Kauerauf 2009, cited in LaRowe et al. 2017
        Cond_sed = 2.5;                                 % W m^-1 K^-1, typical value from LaRowe et al 2017
        waterColumnDepth = 1000;                        % Abritrary water column depth (has very minor impact on solubilities via pressure)
    case 'abyss'
        c = 0.85 * 10^-3;                               % Sediment compaction length scale, from LaRowe et al. 2017
        por_0 = 0.70;                                   % porosity at sediment-water interface (NB, not the very upper cm). Indicative value from Hantschel and Kauerauf 2009, cited in LaRowe et al. 2017
        Cond_sed = 1.7;                                 % W m^-1 K^-1, typical value from LaRowe et al 2017
        waterColumnDepth = 3000;                        % Abritrary water column depth (has very minor impact on solubilities via pressure)
end

a = 3.7980;                                 % This comes from the sensitivity calibration in figure S10 in Tatzel et al., 2022
dz = 0.005;                                 % step size in depth domain in km (strictly valid for first box only; the rest are compressed relative to this).
N = 750;                                        % Total number of nodes in depth domain
dens = 2.6;                                     % density of solid sediment, g/cm3 or t/m3
density_water = 1.030*10^6;                     % g/m3.
[z_real,z,thickness,sed_amt] = defineDepthDomain(dz,N,c,por_0,dens);    % z_real, z, and thickness in m, sed_amt in tons

z_real = z_real(1:end-1) + diff(z_real)/2;      % Midpoint depth of each box, in metres

sedimentDensity = sed_amt./thickness;           % g cm-3, or t/m3
porosity =  1 - sedimentDensity/dens;           % porosity, dimensionless, after accounting for compaction (might later couple this to opalCT or qtz precipitation?), at the box midpoint
SVP = (1-porosity).*thickness;                  % sediment volume profile, m3 of sediment per layer (should be constant, is only here as a check)
porespace = thickness.*porosity;                % what is the actual volume in m3 of available porespace per box?

mass = (porespace * density_water/1000) + sed_amt*1000; % kg of mass in each box. Decreases because the amount of sediment is constant (about 3.3 tons), but the amount of water decreases (from 3.5m3 to 0.5m3).
cumMass = cumsum(mass);

z = [0,z,z(end)];                               % Ghost nodes for numerical schemes
z_real = [0,z_real,z_real(end)];
porespace = [dz*1000,porespace,porespace(end)]; % m3 of porespace per box
thickness = [dz*1000,thickness,thickness(end)]; %
porosity = [1,porosity,porosity(end)];          % extend the length of the vectors by one at either end

Cond_fluid = 0.6;                               % W m^-1 K^-1, from Castelli et al 1974 cited in LaRowe 2017
thermalConductivity = (Cond_sed.^(1-porosity)).*(Cond_fluid.^porosity);     % Geometric mean of sediment and fluid thermal conductivities after Fuchs et al 2013 in LaRowe et al 2017

atmPressure = 1 * 101325;                       % Atmospheric pressure in Pascals (1 = 1 atm)
g = 9.81;                                       % Acceleration due to gravity, 9.81 m/s2

SWIPressure = atmPressure + waterColumnDepth * g * density_water/1000;       % Pressure at the Sediment-Water-Interface, in Pascals. 1000 to convert from g to kg density (kg/m3)
SedPressure = SWIPressure + cumMass * g;        % pressure on water at each box in the sediment column, in Pa.

% ################### Define basic isotope and reaction variables ###################

fractionationModel_O = 'Sharp';                   % can be one of: 'Sharp' [recommended], 'KandE76', 'Meheut2007', or  'BJ73', for Sharp et al., Knauth and Epstein, Meheut et al, or Bottinga and Javoy. See also Pollington et al. 2016

R_SMOW = 2005.20 * 10^-6;                       % ratio of 18/16 in standard mean ocean model, SMOW 0.002
d18_w_i = -1.1;                                 % d18O value of the water, permil
R_bw_O = (d18_w_i/1000 + 1) * R_SMOW;             % 18O/16O ratio of bottom water

% The calibrated isotope abundance ratios of NBS-28 are 0.0507446 for 29Si/28Si and 0.0341465 for 30Si/28Si.

fractionationModel_Si = 'Dupuis';   % can be one of: 'Dupuis' 'Stamm'

beta= (1/28-1/29)/(1/28-1/30); % from Young et al., 2002, alpha29= alpha(30)^beta, equilibrium
R_NBS_3028= 0.0341465;  % ratio of 30/28Si in NBS
R_NBS_2928= 0.0507446;  % ratio of 29/28Si in NBS
d30_w_i= 1;  % d30Si value of the water,permil
R_bw_Si30= (d30_w_i/1000+1) * R_NBS_3028; % 30/28Si of the water
R_bw_Si29= (d30_w_i/1000+1)^ beta * R_NBS_2928; % 29/28Si of the water

R = 8.3144598;                                  % universal gas constant in J/(K*mol)

% SiO2 phase dissolution rate parameters. These are combined to produce a
% rate constant 'k' (=A.exp(-Ea/RT, units = mol/mol kyr-1) and then further
% combined with the equlibrium constants to produce a precipitation rate
% constant. 

[Ea,A] = rateConstantPlots;
A = A * 1000; %convert to units of kyr-1,the pre-exponential factor or frequency factor, which represents the frequency of collisions or attempts at reaction.
Ea = Ea * 1000; % convert to units of J/mol,the activation energy, the minimum energy required for a reaction to occur.

SSA(1) = 75;                                  % m2/mol. (Dixit et al., 2001; Van Cappellen et al., 2002; Van Cappellen and Qiu, 1997a, b; Weiler and Mills, 1965)
SSA(2) = 25;
SSA(3) = 1.38 / 10;

diff_coeff_Si = 0.032 * 1000;    % Diffusion coefficient of dSi in seawater at 25 degrees in m2 kyr-1.  From Wollast and Garrels, updated and corroborated by Rebaneau et al.
% "The first observation is that the diffusion coefficient of dissolved silica measured at salinity 36 and 25 °C in our study, (1.02 ± 0.02) × 10− 5
% cm2 s− 1, is essentially identical to that obtained by Wollast and Garrels (1971), (1.00 ± 0.05) × 10− 5 cm2 s− 1."

sol_dec_coeff = solubilitySuppression(0);           % how much does solibility increase/decrease with detrital mineral content? solubility suppression = ax^2, where x is fraction detral mineral in sample. From Dixit et al data, and others.
% Enter '1' in function to produce a plot

% ################### Define the parameter ensembles to be tested ###################

%%

try
 
    % ################### what should the timestep be?  ###################
    dt = dz * 1000 / v_surf;                  % 500, calculated time step, in kyrs, to move one box.  1000 accounts for unit matching (km * 1000) / m / kyr = kyr
    tmax = max(z)/v_surf;                     % 375250, time to run for, kyr (i.e. one complete sediment profile). All units in m.
    
    SiSed = SVP(1)*dens*10^6*Sed_Si_pc/100/60/dt;        % Amount of sediment in box (m3) x mean density of sediment (t/m3) x fraction SiO2 (t/t) / molar mass (g/mol) x fraction timestep (yr/yr)  = moles of opalA to add to top box per kyr
    
    temp_vec = bwT + z_real*heatFlow./thermalConductivity;          % Temperature in each box, in K
    
%     [a1a, a1b, a1c, a2, a3, a4] = makeProgressPlots(z_real,temp_vec);   % The axes for plotting of progress as we go
    
    % ################### Calculate how fractionation, diffusivity and reaction vary with depth  ###################
    alpha_O = defineFractionation(fractionationModel_O,temp_vec);           % define how fractionation varies as a function of depth (i.e. temperature). Assumed valid for all phases (see Sharp et al.)
    alpha_Si = defineFractionation(fractionationModel_Si,temp_vec); 

    [D_H2O, D_Si] = defineDiffusivities(temp_vec,porosity,diff_coeff_Si);   % define temp- and tortuosity corrected diffusion coefficients for Si and water, units of m2/kyr
    solubilities = defineSiSolubilities(SWIPressure,SedPressure,density_water,temp_vec,sol_dec_coeff,mineralContent,SSA);
    rateConstants = defineRateConstants(A,Ea,R,temp_vec,a,mineralContent/100);

 [massStruct,ratioStruct] = preallocateVariables(z,R_SMOW,d18_w_i,overlying_water_Si_conc,porespace,density_water,solubilities.quartz,R_NBS_2928,R_NBS_3028,d30_w_i,beta); %load(preallocatedVariables)

    nsteps = ceil(tmax/dt);       % 751, how many timesteps requires to simulate the time required to advect the whole sediment column?
    
    % ################### Iterate through time  ###################
    % Here, we increase the number of boxes considered in the mass balance by one every timestep.
    
    i = 0; steadyStateCount = 0; ii = 1;

    
    while i < (nsteps - 1) && steadyStateCount < 20    % Criteria by which to end the simulation       
        
        if  i > 100
            if sum(massStruct.quartz_Si > (0.99 * max(massStruct.opalCT_Si))) > 10  % make sure there's actually more qtz than opalCT
                steadyStateCount = steadyStateCount + 1;
            end
        end
     
        i = i+1;        
        ii = i + 1;
        
        % ################### Reactions ###################
        % 1. Dissolution, maturation and reprecipitation
        
        t_span = [0 dt];
        
        
        C0 = [massStruct.porewater_Si(1:ii)';  
            massStruct.opalA_Si(1:ii)';   
            massStruct.opalCT_Si(1:ii)'; 
            massStruct.quartz_Si(1:ii)';  
            massStruct.porewater_16O(1:ii)'; 
            massStruct.porewater_18O(1:ii)';
            massStruct.opalA_16O(1:ii)'; 
            massStruct.opalA_18O(1:ii)'; 
            massStruct.opalCT_16O(1:ii)';
            massStruct.opalCT_18O(1:ii)';
            massStruct.quartz_16O(1:ii)'; 
            massStruct.quartz_18O(1:ii)' ;
            massStruct.porewater_28Si(1:ii)';
            massStruct.porewater_30Si(1:ii)'; 
            massStruct.opalA_28Si(1:ii)'; 
            massStruct.opalA_30Si(1:ii)'; 
            massStruct.opalCT_28Si(1:ii)'; 
            massStruct.opalCT_30Si(1:ii)';
            massStruct.quartz_28Si(1:ii)';
            massStruct.quartz_30Si(1:ii)'
           ];
    
        % the fundamental mass-balance step. Includes dissolution and
        % precipitation, advection (of all phases) and diffusion (of
        % solutes).
        
        tic
        [TOUT,COUT] = ode15s(@(t,C) ...
            SiO2ReactionAdvectionODEs(C,solubilities,rateConstants,ii,v_surf,dz,SiSed,D_Si,D_H2O,alpha_O,porespace, alpha_Si, beta,fractionation_si),...
            t_span,C0);
        
        tt = toc; disp(strcat("Step = ",num2str(i),{'; time taken = '},num2str(tt),{' seconds'}))
        
        
        % unpack the results:
        n = size(COUT,2) / 20;  % the width of the matrix, divided by the number of equations (12x silicon, 8 x oxygen)
        
        newPorewater_Si = COUT(end,1:n);
        newOpalA_Si = COUT(end,n*1+1:n*2);
        newOpalCT_Si = COUT(end,n*2+1:n*3);
        newQuartz_Si = COUT(end,n*3+1:n*4);
        
        newPorewater_16O = COUT(end,n*4+1:n*5);
        newPorewater_18O = COUT(end,n*5+1:n*6);
        newOpalA_16O = COUT(end,n*6+1:n*7);
        newOpalA_18O = COUT(end,n*7+1:n*8);
        newOpalCT_16O = COUT(end,n*8+1:n*9);
        newOpalCT_18O = COUT(end,n*9+1:n*10);
        newQuartz_16O = COUT(end,n*10+1:n*11);
        newQuartz_18O = COUT(end,n*11+1:n*12);

        newPorewater_28Si = COUT(end,n*12+1:n*13);
        newPorewater_30Si = COUT(end,n*13+1:n*14);
        newOpalA_28Si = COUT(end,n*14+1:n*15);
        newOpalA_30Si = COUT(end,n*15+1:n*16);
        newOpalCT_28Si = COUT(end,n*16+1:n*17);
        newOpalCT_30Si = COUT(end,n*17+1:n*18);
        newQuartz_28Si = COUT(end,n*18+1:n*19);
        newQuartz_30Si = COUT(end,n*19+1:n*20);

        % update variables
        massStruct.porewater_Si(1:ii) = newPorewater_Si;
        massStruct.opalA_Si(1:ii) = newOpalA_Si;
        massStruct.opalCT_Si(1:ii) = newOpalCT_Si;
        massStruct.quartz_Si(1:ii) = newQuartz_Si;
        
        massStruct.porewater_16O(1:ii) = newPorewater_16O;
        massStruct.porewater_18O(1:ii) = newPorewater_18O;
        massStruct.opalA_16O(1:ii) = newOpalA_16O;
        massStruct.opalA_18O(1:ii) = newOpalA_18O;
        massStruct.opalCT_16O(1:ii) = newOpalCT_16O;
        massStruct.opalCT_18O(1:ii) = newOpalCT_18O;
        massStruct.quartz_16O(1:ii) = newQuartz_16O;
        massStruct.quartz_18O(1:ii) = newQuartz_18O;

        massStruct.porewater_28Si(1:ii) = newPorewater_28Si;
        massStruct.porewater_30Si(1:ii) = newPorewater_30Si;
        massStruct.opalA_28Si(1:ii) = newOpalA_28Si;
        massStruct.opalA_30Si(1:ii) = newOpalA_30Si;
        massStruct.opalCT_28Si(1:ii) = newOpalCT_28Si;
        massStruct.opalCT_30Si(1:ii) = newOpalCT_30Si;
        massStruct.quartz_28Si(1:ii) = newQuartz_28Si;
        massStruct.quartz_30Si(1:ii) = newQuartz_30Si;
        
       
        porewater_Si_conc = massStruct.porewater_Si./(porespace*1000);
        porewater_Si_conc(ii+1) = porewater_Si_conc(ii);
        massStruct.porewater_Si(ii+1) = porewater_Si_conc(ii+1) * (porespace(ii+1)*1000);
        
        % If masses of opalA or opalCT get too low, set to zero, to avoid odd numerical problems:
        massStruct.opalA_Si([massStruct.opalA_Si] < 1e-10) = 0;
        massStruct.opalA_16O([massStruct.opalA_Si] < 1e-10) = 0;
        massStruct.opalA_18O([massStruct.opalA_Si] < 1e-10) = 0;
         massStruct.opalA_28Si([massStruct.opalA_Si] < 1e-10) = 0;
         massStruct.opalA_30Si([massStruct.opalA_Si] < 1e-10) = 0;
  
        massStruct.opalCT_Si([massStruct.opalCT_Si] < 1e-10) = 0;
        massStruct.opalCT_16O([massStruct.opalCT_Si] < 1e-10) = 0;
        massStruct.opalCT_18O([massStruct.opalCT_Si] < 1e-10) = 0;
         massStruct.opalCT_28Si([massStruct.opalCT_Si] < 1e-10) = 0;
         massStruct.opalCT_30Si([massStruct.opalCT_Si] < 1e-10) = 0;


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

%         if rem(i,50) == 0 % Figure plotting code section
%             updateProgressPlots(z_real,massStruct,ratioStruct,solubilities,dt,porespace,a1a, a1b, a1c, a2, a3,a4, ii, SVP(1)*dens*10^6)
%         end

    end
    
    if isfolder('newresults') == 1
        cd newresults
    else 
        mkdir('newresults');
        cd newresults
    end
    
%     close all
    
    filename = strcat('ChertModelOutput',datestr(datetime('now')));
    save(filename)
    
    cd ..
    
    [results] = findResults(massStruct,ratioStruct,z_real,temp_vec);
    disp(['Depth when quartz is >99.5% of all SiO2 = ',num2str(round(results(1),1)),'m; chert d18O = ',num2str(round(results(2),2)),' permil; chert d30Si = ',num2str(round(results(4),4)),' permil;  temperature = ',num2str(round(results(3)-273.15,1)),' deg C'])

catch ME
    [~,errorInFile,~] = fileparts(ME.stack(1).file);
    warning(['some issue with iteration on line ',num2str(ME.stack(1).line),' of function ',errorInFile])
end


disp('Finished!')
