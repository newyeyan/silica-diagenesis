function [D_H2O, D_Si] = defineDiffusivities(temp_vec,porosity,diff_coeff_Si)
% PJF 2022
% calculate water and dissolved Si diffusion coefficients for the sediment
% column (as a function of temp), implementing a Boudreau 1996 style
% correction for tortousity. Units of m2/kyr throughout. 

   
    D_m = -925.1679617;                                                         % Gradient of the temperature/water self-diffusivity coefficient relationship, from Tanaka (1975)
    D_b = -5.55858328;                                                          % Intercept of the temperature/water self-diffusivity coefficient relationship, from Tanaka (1975)
    D_H2O = 10.^((D_m.* 1./temp_vec)+D_b);                                      % diffusivity of water as a function of T, derived from Tanaka 1974; see also Easteal et al. (1984) for potential isotope fractionation during diffusion
    tortuosity = 1 - (log(porosity.^2)); tortuosity(1) = tortuosity(2);         % Tortuosity squared. Workaround for point 1 reuired because otherwise diffusion too large at sediment water interface
    D_H2O = (D_H2O .* porosity) ./ (tortuosity);                                % effective diffusivity, after correcting for tortuosity after Boudreau 1996
    D_H2O = D_H2O * 60 * 60 * 24 * 365 * 1000;                                  % Convert to units of m2/kyr.  NB: when plotted, has an odd shape because of competing effects of porosity/tortuosity and T.

  % Si diffusivity as a function of T and viscosity
    seawater_viscosity =  0.1899 ./ (1 - 0.01203 *temp_vec + (3.307*10^-5)*temp_vec.^2) * 0.001787; % modified from the eqn of Hardie to get into units of T in K
    V_at_298 = 0.1899 ./ (1 - 0.01203 *298 + (3.307*10^-5)*298^2) * 0.001787;                       % Seawater viscosity at 298 K
    D_uncorr = (diff_coeff_Si  * V_at_298 / 298) .* temp_vec ./ seawater_viscosity;                 % Einstein-Stokes equation correction, for T and viscoscity. See e.g. Rebreanu et al. 2008 Mar Chem for details.
    D_Si = (porosity.*D_uncorr ) ./ tortuosity;                                                     % Boudreau correction for tortuosity. Empirical, but other formulations (Ullman and Aller/Ridgwell) based on formation factor) give similar results
    
end