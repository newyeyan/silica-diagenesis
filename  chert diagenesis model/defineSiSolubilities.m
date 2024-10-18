function [solubilities] = defineSiSolubilities(SWIPressure,SedPressure,density_water,temp_vec,sol_dec_coeff,mineralContent,SSA)

% PJF 2021: Define the solubilities of opal-A, opal-CT and quartz as a
% function of temperautre and pressure. Use expression form von Damm et al
% for quartz, and the expression given in Sarmiento and Gruber for opal-A
% (this is itself a combination of multiple studies). No reliable data
% exists for opal-CT solubility in seawater, so take a point lying 25% of
% the way between quartz and opal A. This gives a reasoable approximation
% of the opal-CT solubility in the freshwater compilation of Gunnarsson and
% Arnursson. Note that opal-A solubility includes the effect of decreasing
% solubility given by the presence of detrital aluminosilciates, and that
% this is therefore included (albeit to a lesser extent) in the opal-CT
% solubility estimate. 

% These all give solubilities in units of moles/litre.

% von Damm expression for quartz
SedPressure = [SWIPressure,SedPressure,SedPressure(end)]/100000;    % Pore pressure, in bars, to agree with units in von Damm et al. 1991 AJS.
vD_a = -2.32888;                                                    % Coefficients for best fit to data from von Damm et al. 1991
vD_b = 1.79547;
vD_c = -2263.62;
vD_d = 0.00407350;
vD_e = 0.0398808;

SolQtz = exp(vD_a + vD_b*log(density_water/1000000) + (vD_c + vD_d .* temp_vec.^2) .* temp_vec.^-1 + (vD_e * SedPressure .* temp_vec.^-1)); % from von Damm et al. the 'a' term contains the impact of NaCl...

% Dixit et al/Sarmiento and Gruber expression for opalA/Biogenic silica
SolASi = 2.754 * exp(1./temp_vec .* (-2229 ...
    - 3.688 .* 10^-3 .* 1500 + 3.688 .* 10^-3 .* SSA(1) ...
    + 0.2 .* SedPressure./1.013 - 2.7*10^-4 .* (SedPressure/1.013).^2 + 1.46*10^-7 .* (SedPressure./1.013).^3));



SolASi = SolASi - SolASi * sol_dec_coeff * (mineralContent/100)^2; % detrital mineral influence.

SolCT = SolQtz + (SolASi - SolQtz)/4;

solubilities.opalA = SolASi;
solubilities.opalCT = SolCT;
solubilities.quartz = SolQtz;






%% plot solubilities
% 
% fig2 = figure;
% fig2.Units = 'centimeters';
% fig2.PaperOrientation = 'landscape';
% fig2.Color = [1 1 1];
% 
% fig2.Position = [0 0 13 13];
% 
% % 1a = 795 concentrations (porewater, solid phases)
% axes1a = axes;
% axes1a.Units = 'centimeters';
% axes1a.Position = [3 3 8 8];
% axes1a.Box = 'on';
% axes1a.XColor = 'k';
% axes1a.YColor = 'k';
% axes1a.Color = 'none';
% axes1a.XLabel.String = 'Temperature (K)';
% axes1a.YLabel.String = 'Concentration (\muM)';
% axes1a.YLim = [0 2000];
% hold on
% 
% p1= plot(temp_vec,solubilities.opalA*10^6,'r-');
% p2= plot(temp_vec,solubilities.opalCT*10^6,'b--');
% p3= plot(temp_vec,solubilities.quartz*10^6,'k-');
% 
% 
% l = legend([p1,p2,p3],{'OpalA','OpalCT','Quartz'});
% l.Location = 'southeast';
% l.Color = 'none';
% l.Box = 'off';

end