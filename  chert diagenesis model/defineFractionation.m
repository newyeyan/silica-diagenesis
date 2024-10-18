function alpha = defineFractionation(fractionationModel,temp_vec)

% PJF 2021: A selection of different SiO2/H2O oxygen isotope/temperature
% calibrations from the literature.

switch fractionationModel
    case 'Meheut2007'
        x = (10^6)./temp_vec.^2;
        fractionation = [-9.6341 + 4.848*x(temp_vec < 403.15) - (0.0382 * (x(temp_vec < 403.15)).^2) + (0.000397* (x(temp_vec < 403.15)).^3),...
            -2.9548 + (1.342*x(temp_vec >= 403.15)) + (0.6062 * (x(temp_vec >= 403.15)).^2 )+ (-0.040638* (x(temp_vec >= 403.15)).^3)];
        
    case 'KandE76'
        fractionation = 3.09 * 10^6 .* temp_vec.^(-2) -3.29;    % the empirical relationship of Knauth and Epstein.
        
    case 'BJ73'
        A = -3.7; B = 4.1 * 10^6;
        fractionation = A + B./(temp_vec).^2;
        
    case 'Sharp'
        A = 4.28*10^6; B = 3.5*10^3;
        fractionation = A./(temp_vec).^2 - B./temp_vec;

    case 'Dupuis'  %  get from the linear relation
        A = 0.267*10^6; B = -0.9;
        fractionation = A./(temp_vec).^2 + B;
        
    case 'Stamm' %  get from the linear relation
        A = 0.133*10^6; B = -1;
        fractionation = A./(temp_vec).^2 + B;
        
    otherwise
        disp('unknown T-fractionation model - ending'), return
end

alpha = fractionation/1000 + 1;                                             % the fractionation expressed as a factor (alpha notation)