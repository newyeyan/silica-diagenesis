function [a1a, a1b, a1c, a2, a3, a4] = makeProgressPlots(z_real,temp_vec)

% PJF September 2021. Just generates a few axes for plotting later. First
% shows concentrations of dSi in porewater along with the solubility lines
% of the three SiO2 phases. Second shows concentrations, in wt%, of the
% three SiO2 phases, and third shows the d81O values of SiO2 and porewater
% dSi. These are very rough, and just for monitoring if everything is
% working more-or-less as planned...
 
    
    fig1 = figure('Name','Porewater concentration, Solubility, Temperature'); 
    a1a = axes; fig1.Position = [1   535   560   420];
    a1b = axes; a1b.Position = a1a.Position; a1b.Color = 'none';
    a1b.YAxisLocation = 'right'; a1b.YLabel.String = 'Temperature (deg C)'; a1b.XTick = [];   hold on
    a1c = axes; a1c.Position = a1a.Position;
    
    fig2 = figure('Name','Solid phase concentrations'); 
    a2 = axes; 
    fig2.Position = [562   535   560   420];

    fig3 = figure('Name','Oxygen isotope ratios'); 
    a3 = axes; 
    fig3.Position = [1121         535         560         420];

    fig4 = figure('Name','Silicon isotope ratios'); 
    a4 = axes; 
    fig4.Position = [1121         1         560         420];

    plot(a1b,z_real,temp_vec-273.15,'r.')
    axes(a1b),   tempLabel = text((a1b.XLim(2) - a1b.XLim(1))/2, temp_vec(ceil(length(temp_vec)/2))-273.15,'Temperature');
    tempLabel.Rotation = 40; tempLabel.Color = 'r'; tempLabel.FontWeight = 'bold'; tempLabel.HorizontalAlignment = 'right';
