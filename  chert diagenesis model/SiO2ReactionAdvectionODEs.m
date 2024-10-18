function dMdt = SiO2ReactionAdvectionODEs(M,solubilities,rateConstants,j,v,dz,SiSed,D_Si,D_H20,alpha_O,porespace, alpha_Si,beta,fractionation_si)

% revised from PJF 2021. This is where the O and Si mass balances are calculated.
% Essentially, this function is defining the time rate of change for each
% phase in each box as dM/dt = reaction + advection + diffusion. Serves as
% input to Matlab ODE solver ode15s, (designed for stiff systems). Can get
% slow when number of boxes increases, beacuse the number of equations to
% solve is n*20, where n is the number of sediment boxes.

% M = masses that we're solving the equations for. A vector of length 20*j distributed as: 
%     M(1) = porewater Si,  M(0*j+1:j*1)
%     M(2) = opalA Si,      M(1*j+1:j*2)
%     M(3) = opalCT Si,     M(2*j+1:j*3)
%     M(4) = quartz Si,     M(3*j+1:j*4)
%     M(5) = porewater 16O, M(4*j+1:j*5)
%     M(6) = porewater 18O, M(5*j+1:j*6)
%     M(7) = opalA 16O,     M(6*j+1:j*7)
%     M(8) = opalA 18O,     M(7*j+1:j*8)
%     M(9) = opalCT 16O,    M(8*j+1:j*9)
%     M(10) = opalCT 18O,   M(9*j+1:j*10)
%     M(11) = quartz 16O,   M(10*j+1:j*11)
%     M(12) = quartz 18O,   M(11*j+1:j*12)
%     M(13) = porewater 28Si,  M(12*j+1:j*13)
%     M(14) = porewater 30Si,  M(13*j+1:j*14)
%     M(15) = opalA 28Si,      M(14*j+1:j*15)
%     M(16) = opalA 30Si,      M(15*j+1:j*16)
%     M(17) = opalCT 28Si,     M(16*j+1:j*17)
%     M(18) = opalCT 30Si,     M(17*j+1:j*18)
%     M(19) = quartz 28Si,     M(18*j+1:j*19)
%     M(20) = quartz 30Si,     M(19*j+1:j*20)

% j = timestep; total number of boxes considered (one timestep is the time
% required to advect sediment the length of one box.


M(M<0) = 0; % in case of odd numerical issues
lowOpalA = M(1*j+2:j*2) < 1e-12 & M(1*j+2:j*2) > 0; % total Si in opal A
M(1*j+1+find(lowOpalA)) = 0;    % set opal-A to 0 if it's less than 10^-12 (i.e. practically zero).
M(6*j+1+find(lowOpalA)) = 0;    % and the corresponding oxygen isotopes too...
M(7*j+1+find(lowOpalA)) = 0;
M(14*j+1+find(lowOpalA)) = 0;
M(15*j+1+find(lowOpalA)) = 0;

dz = dz * 1000; % km to m conversion,
s_Si = (D_Si)'/(dz^2);  % s terms in FTCS diffusion scheme. Unit check: dz = m, D_Si = m2/kyr (tortuosity etc corrected), and all time units are in kyr.
s_H2O = (D_H20)'/(dz^2);

% Conc represents the mass of phase j at timestep n and spatial node i
Conc_Si_expanded = [M(1:j);M(j)]; %./ [porespace(1:j)';porespace(j)];  add extra nodes for the porewaters (at top = ocean water; at bottom = no-flow BC, so lower-most value).
waterOxygen_expanded = [M(4*j+1:j*5) + M(5*j+1:j*6); M(j*5) + M(j*6)];
Conc_16O_expanded = [M(4*j+1:j*5);M(j*5)]; % ./ waterOxygen_expanded;
Conc_18O_expanded = [M(5*j+1:j*6);M(j*6)]; % ./ waterOxygen_expanded;

Conc_28Si_expanded = [M(12*j+1:j*13);M(j*13)] ;
Conc_30Si_expanded = [M(13*j+1:j*14);M(j*14)] ;

% Conc_16O_expanded = [M(4*j+1:j*5);M(j*5)] ./ [porespace(1:j)';porespace(j)];
% Conc_18O_expanded = [M(5*j+1:j*6);M(j*6)] ./ [porespace(1:j)';porespace(j)];

% Reaction (diss/precip) Ultimately derived from eqn of the form given in e.g. Kamatani et al 1980. 
dOpalAdt = - M(j+1:2*j) ./ porespace(1:j)' .* (rateConstants.k1(1:j)' ./ solubilities.opalA(1:j)') .* (solubilities.opalA(1:j)' - M(1:j)./(porespace(1:j)'*1000)); 
dOpalCTdt = - M(2*j+1:j*3) ./ porespace(1:j)' .* (rateConstants.k2(1:j)' ./ solubilities.opalCT(1:j)') .* (solubilities.opalCT(1:j)' - M(1:j)./(porespace(1:j)'*1000));
dQuartzdt = - M(3*j+1:j*4) ./ porespace(1:j)' .* (rateConstants.k3(1:j)' ./ solubilities.quartz(1:j)') .* (solubilities.quartz(1:j)' - M(1:j)./(porespace(1:j)'*1000));

dOpalAdt((dOpalAdt + M(j+1:j*2)) < 0) = - M(j+find((dOpalAdt + M(j+1:j*2)) < 0)); 
dOpalCTdt((dOpalCTdt + M(2*j+1:j*3)) < 0) = - M(2*j+find((dOpalCTdt + M(2*j+1:j*3)) < 0));
dQuartzdt((dQuartzdt + M(3*j+1:j*4)) < 0) = - M(3*j+find((dQuartzdt + M(3*j+1:j*4)) < 0));

concs_Si = Conc_Si_expanded(2:j); concs_Si_minus1 = Conc_Si_expanded(1:j-1); concs_Si_plus1 = Conc_Si_expanded(3:end);
concs_16O = Conc_16O_expanded(2:j); concs_16O_minus1 = Conc_16O_expanded(1:j-1); concs_16O_plus1 = Conc_16O_expanded(3:end);
concs_18O = Conc_18O_expanded(2:j); concs_18O_minus1 = Conc_18O_expanded(1:j-1); concs_18O_plus1 = Conc_18O_expanded(3:end);
concs_28Si = Conc_28Si_expanded(2:j); concs_28Si_minus1 = Conc_28Si_expanded(1:j-1); concs_28Si_plus1 = Conc_28Si_expanded(3:end);
concs_30Si = Conc_30Si_expanded(2:j); concs_30Si_minus1 = Conc_30Si_expanded(1:j-1); concs_30Si_plus1 = Conc_30Si_expanded(3:end);

% Porewater total Si
dMdt1(1:j,1) = - dOpalAdt - dOpalCTdt - dQuartzdt;

% Conc represents the mass of phase j at timestep n and spatial node i
% c(n+1,j,i) = s * c(n,j,i-1) + (1-2s) * c(n,j,i) + s* c(n,j,i+1),
% = c(n,j,i) + s * (c(n,j,i-1) - 2 * c(n,j,i) + c(n,j,i+1)).
% diffusion of DSi (FTCS)
dMdt1(2:j,1) = dMdt1(2:j,1) + s_Si(2:j).*(concs_Si_minus1 - 2*concs_Si + concs_Si_plus1);  

% opalA total Si
dMdt2(1:j,1) = dOpalAdt;   % diss or precip of opalA, if negative, dissolve, with the value of opal A, 
% if positive, precipitate with fractionation from porewater

% opalCT total Si
dMdt3(1:j,1) = dOpalCTdt;   % diss or precip of opalCT

% quartz total Si
dMdt4(1:j,1) = dQuartzdt; % diss or precip of Qtz

% opalA oxygen, if negative, dissolve, with the value of opal A, if positive, precipitate with fractionation
% opalA O16
dMdt7(1:j,1) = (dOpalAdt .* (dOpalAdt < 0) ./ (1 + M(7*j+1:8*j)./M(6*j+1:7*j))...
    + dOpalAdt .* (dOpalAdt > 0) ./ (1 + alpha_O(1:j)' .*  M(5*j+1:6*j)./M(4*j+1:5*j)) ) * 2;
dMdt7(isnan(dMdt7)) = 0;
% opalA O18
dMdt8(1:j,1) = dOpalAdt*2 - dMdt7;

% opalCT oxygen
% opalCT O16
dMdt9(1:j,1) = (dOpalCTdt .* (dOpalCTdt < 0) ./ (1 + M(9*j+1:10*j)./M(8*j+1:9*j))...
    + dOpalCTdt .* (dOpalCTdt > 0) ./ (1 + alpha_O(1:j)' .*  M(5*j+1:6*j)./M(4*j+1:5*j)) ) * 2;
dMdt9(isnan(dMdt9)) = 0;
% opalCT O18
dMdt10(1:j,1) = dOpalCTdt * 2 - dMdt9;

if sum(M(8*j+1:9*j) < 0) > 0
    j
end

% quartz oxygen O16
dMdt11(1:j,1) = (dQuartzdt .* (dQuartzdt < 0) ./ (1 + M(11*j+1:12*j)./M(10*j+1:11*j))...
    + dQuartzdt .* (dQuartzdt > 0) ./ (1 + alpha_O(1:j)' .*  M(5*j+1:6*j)./M(4*j+1:5*j)) ) * 2;
dMdt11(isnan(dMdt11)) = 0;
% quartz oxygen O18
dMdt12(1:j,1) = dQuartzdt*2 - dMdt11;

if sum(M(10*j+1:11*j) < 0) > 0
    j
end

% porewater oxygen
dMdt5 = - dMdt7 - dMdt9 - dMdt11;
dMdt6 = - dMdt8 - dMdt10 - dMdt12;

% diffusion of porewater
dMdt5(2:j,1) = dMdt5(2:j,1) + s_H2O(2:j).*(concs_16O_minus1 - 2*concs_16O + concs_16O_plus1); %.*waterOxygen_expanded(2:j);
dMdt6(2:j,1) = dMdt6(2:j,1) + s_H2O(2:j).*(concs_18O_minus1 - 2*concs_18O + concs_18O_plus1); %.*waterOxygen_expanded(2:j);

%
% opal A 28Si
dMdt15(1:j,1) = dOpalAdt .* (dOpalAdt < 0) .* M(14*j+1:15*j) ./M(1*j+1:j*2) ...
    + dOpalAdt .* (dOpalAdt > 0) ./ (1 + alpha_Si(1:j)' .*  M(13*j+1:14*j)./M(12*j+1:13*j) ...
    + alpha_Si(1:j)' .^beta .*  (M(0*j+1:j*1)-M(13*j+1:14*j)- M(12*j+1:13*j))./M(12*j+1:13*j)); 
dMdt15(isnan(dMdt15)) = 0;
% % opal A 30Si
dMdt16(1:j,1) = dOpalAdt .* (dOpalAdt < 0) .* M(15*j+1:16*j) ./M(1*j+1:j*2) ...
    + dMdt15(1:j,1).* (dOpalAdt > 0) .*  alpha_Si(1:j)' .*  M(13*j+1:14*j)./M(12*j+1:13*j);  %opal A 30Si
dMdt16(isnan(dMdt16)) = 0;
% 
% % opalCT 28Si
dMdt17(1:j,1) = dOpalCTdt .* (dOpalCTdt < 0) .* M(16*j+1:17*j) ./M(2*j+1:j*3) ...
    + dOpalCTdt .* (dOpalCTdt > 0) ./ (1 + alpha_Si(1:j)' .*  M(13*j+1:14*j)./M(12*j+1:13*j) ...
    + alpha_Si(1:j)' .^beta .*  (M(0*j+1:j*1)- M(13*j+1:14*j)- M(12*j+1:13*j))./M(12*j+1:13*j)); 
dMdt17(isnan(dMdt17)) = 0;
% %opal CT 30Si
dMdt18(1:j,1) = dOpalCTdt .* (dOpalCTdt < 0) .* M(17*j+1:18*j) ./M(2*j+1:j*3) ...
    + dMdt17(1:j,1).* (dOpalCTdt > 0) .*  alpha_Si(1:j)' .*  M(13*j+1:14*j)./M(12*j+1:13*j);  %opal CT 30Si
% 
dMdt18(isnan(dMdt18)) = 0;
if sum(M(16*j+1:17*j) < 0) > 0
    j
end
% 
% 
% % quartz 28Si
dMdt19(1:j,1) = dQuartzdt .* (dQuartzdt < 0) .* M(18*j+1:19*j) ./M(3*j+1:j*4) ...
    + dQuartzdt .* (dQuartzdt > 0) ./ (1 + alpha_Si(1:j)' .*  M(13*j+1:14*j)./M(12*j+1:13*j) ...
    + alpha_Si(1:j)' .^beta .*  (M(0*j+1:j*1)-M(13*j+1:14*j)- M(12*j+1:13*j))./M(12*j+1:13*j)); 
dMdt19(isnan(dMdt19)) = 0;
% % quartz 30Si
dMdt20(1:j,1) = dQuartzdt .* (dQuartzdt < 0) .* M(19*j+1:20*j) ./M(3*j+1:j*4) ...
    + dMdt19(1:j,1).* (dQuartzdt > 0) .*  alpha_Si(1:j)' .*  M(13*j+1:14*j)./M(12*j+1:13*j);  
dMdt20(isnan(dMdt20)) = 0;
% 
if sum(M(18*j+1:19*j) < 0) > 0
    j
end
% 
% % porewater silicon
dMdt13 = - dMdt15 - dMdt17 - dMdt19; %28Si
dMdt14 = - dMdt16 - dMdt18 - dMdt20; %30Si

% 
% % diffusion of porewater 30Si and 28Si
dMdt13(2:j,1) = dMdt13(2:j,1) + s_Si(2:j).*(concs_28Si_minus1 - 2*concs_28Si + concs_28Si_plus1); %.*porespace(2:j)';
dMdt14(2:j,1) = dMdt14(2:j,1) + s_Si(2:j).*(concs_30Si_minus1 - 2*concs_30Si + concs_30Si_plus1); %.*porespace(2:j)';

% advection
% dMdt1(2:j) = dMdt1(2:j) + v/dz * (M((1-1)*j+1:1*j-1) - M((1-1)*j+2:1*j));
dMdt2(2:j) = dMdt2(2:j) + v/dz * (M((2-1)*j+1:2*j-1) - M((2-1)*j+2:2*j));
dMdt3(2:j) = dMdt3(2:j) + v/dz * (M((3-1)*j+1:3*j-1) - M((3-1)*j+2:3*j));
dMdt4(2:j) = dMdt4(2:j) + v/dz * (M((4-1)*j+1:4*j-1) - M((4-1)*j+2:4*j));
% dMdt5(2:j) = dMdt5(2:j) + v/dz * (M((5-1)*j+1:5*j-1) - M((5-1)*j+2:5*j));
% dMdt6(2:j) = dMdt6(2:j) + v/dz * (M((6-1)*j+1:6*j-1) - M((6-1)*j+2:6*j));
dMdt7(2:j) = dMdt7(2:j) + v/dz * (M((7-1)*j+1:7*j-1) - M((7-1)*j+2:7*j));
dMdt8(2:j) = dMdt8(2:j) + v/dz * (M((8-1)*j+1:8*j-1) - M((8-1)*j+2:8*j));
dMdt9(2:j) = dMdt9(2:j) + v/dz * (M((9-1)*j+1:9*j-1) - M((9-1)*j+2:9*j));
dMdt10(2:j) = dMdt10(2:j) + v/dz * (M((10-1)*j+1:10*j-1) - M((10-1)*j+2:10*j));
dMdt11(2:j) = dMdt11(2:j) + v/dz * (M((11-1)*j+1:11*j-1) - M((11-1)*j+2:11*j));
dMdt12(2:j) = dMdt12(2:j) + v/dz * (M((12-1)*j+1:12*j-1) - M((12-1)*j+2:12*j));
dMdt13(2:j) = dMdt13(2:j) + v/dz * (M((13-1)*j+1:13*j-1) - M((13-1)*j+2:13*j));
dMdt14(2:j) = dMdt14(2:j) + v/dz * (M((14-1)*j+1:14*j-1) - M((14-1)*j+2:14*j));
dMdt15(2:j) = dMdt15(2:j) + v/dz * (M((15-1)*j+1:15*j-1) - M((15-1)*j+2:15*j));
dMdt16(2:j) = dMdt16(2:j) + v/dz * (M((16-1)*j+1:16*j-1) - M((16-1)*j+2:16*j));
dMdt17(2:j) = dMdt17(2:j) + v/dz * (M((17-1)*j+1:17*j-1) - M((17-1)*j+2:17*j));
dMdt18(2:j) = dMdt18(2:j) + v/dz * (M((18-1)*j+1:18*j-1) - M((18-1)*j+2:18*j));
dMdt19(2:j) = dMdt19(2:j) + v/dz * (M((19-1)*j+1:19*j-1) - M((19-1)*j+2:19*j));
dMdt20(2:j) = dMdt20(2:j) + v/dz * (M((20-1)*j+1:20*j-1) - M((20-1)*j+2:20*j));

% Boundary conditions: no reaction in box1 (this is not "sediment"...)
dMdt1(1) = 0;
dMdt2(1) = 0;
dMdt3(1) = 0;
dMdt4(1) = 0;
dMdt5(1) = 0;
dMdt6(1) = 0;
dMdt7(1) = 0;
dMdt8(1) = 0;
dMdt9(1) = 0;
dMdt10(1) = 0;
dMdt11(1) = 0;
dMdt12(1) = 0;
dMdt13(1) = 0;
dMdt14(1) = 0;
dMdt15(1) = 0;
dMdt16(1) = 0;
dMdt17(1) = 0;
dMdt18(1) = 0;
dMdt19(1) = 0;
dMdt20(1) = 0;

% Boundary conditions continued: sedimentation of opalA onto first box of sediment:
dMdt2(2) = dMdt2(2) + SiSed; % opal A, si, Amount of sediment in box (m3) x mean density of sediment (t/m3) x fraction SiO2 (t/t) / molar mass (g/mol) x fraction timestep (yr/yr)  = moles of opalA to add to top box per kyr
dMdt7(2) = dMdt7(2) + (SiSed / (1 + alpha_O(1) *  M(5*j+1)./M(4*j+1))) * 2; % 16O
dMdt8(2) = dMdt8(2) + SiSed*2 - ( (SiSed / (1 + alpha_O(1) *  M(5*j+1)./M(4*j+1))) * 2); % 18O
% 
 alpha_Si(1) = fractionation_si ; % initial opal A 
dMdt15(2) = dMdt15(2) + SiSed / (1 + alpha_Si(1) *  M(13*j+1)./M(12*j+1) + alpha_Si(1) .^beta *  (M(0*j+1)-M(13*j+1)-M(12*j+1))./M(12*j+1) );
dMdt16(2) = dMdt16(2) + SiSed .* alpha_Si(1) *  M(13*j+1)./M(12*j+1) / (1 + alpha_Si(1) *  M(13*j+1)./M(12*j+1) + alpha_Si(1) .^beta *  (M(0*j+1)-M(13*j+1)-M(12*j+1))./M(12*j+1)) ;

dMdt = [dMdt1;dMdt2;dMdt3;dMdt4;dMdt5;dMdt6;dMdt7;dMdt8;dMdt9;dMdt10;dMdt11;dMdt12;dMdt13;dMdt14;dMdt15;dMdt16;dMdt17;dMdt18;dMdt19;dMdt20];
end