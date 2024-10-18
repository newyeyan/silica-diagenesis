function [Ea,A] = rateConstantPlots

% PJF 2021: derive rate constants for dissolution of opal-A and quartz in
% seawater as a fucntion of temperature. We're trying to get the
% preexponential term and the activation energy of an Arrhenius expression.
% opalA data from Icenhower and Dover 2000, quartz from Dove and Crerar
% 1990, in both cases with the dataset filtered to include only solutions
% with ionic strenght similar to seawater. Beacuse no good experimental
% data exists for opal-CT, simply take it as the mean of opal-A and quartz.

R = 8.3144598; % j K-1 mol-1, universal gas constant.

%% OpalA
cd auxFiles/inputData

ID2000 = load('IcenhowerDove2000.txt');

opalA_temps = ID2000(:,1);
opalA_rates = ID2000(:,2);
opalA_rates = opalA_rates * 60 * 60 * 24 * 365;

opalA_rates = log(opalA_rates); % express as natural log rate
opalA_temps = opalA_temps + 273.15; % convert to Kelvin
opalA_temps = 1./opalA_temps;

[xData, yData] = prepareCurveData( opalA_temps, opalA_rates );
ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft );
xVals = linspace(min(opalA_temps)*0.75,max(opalA_temps)*1.25,50);
CI = confint(fitresult);
PI = predint(fitresult,xVals,0.95);


Ea(1) = -fitresult.p1*R/1000;
A(1) = exp(fitresult.p2);
errOnEa(1) = -CI(1,1) * R / 1000 - Ea(1);


%% Quartz

DC1990 = load('DoveCrerar1990.txt');

quartz_temps = DC1990(:,1);
quartz_rates = DC1990(:,2);


quartz_rates = quartz_rates * 60 * 60 * 24 * 365;
quartz_rates = log(quartz_rates); % express as natural log rate
quartz_temps = quartz_temps + 273.15; % convert to Kelvin
quartz_temps = 1./quartz_temps;

[xData, yData] = prepareCurveData( quartz_temps, quartz_rates );
ft = fittype( 'poly1' );
[fitresult2, gof2] = fit( xData, yData, ft );
xVals2 = linspace(min(quartz_temps)*0.75,max(quartz_temps)*1.25,50);
CI2 = confint(fitresult2);
PI2 = predint(fitresult2,xVals2,0.95);


Ea(3) = -fitresult2.p1*R/1000;
A(3) =  exp(fitresult2.p2);
errOnEa(3) = -CI2(1,1) * R / 1000 - Ea(3);


%% opalCT
Ea(2) = mean([Ea(1),Ea(3)]);
A(2) = mean([A(1),A(3)]);


% %% plotting
% 
% fig2 = figure;
% fig2.Units = 'centimeters';
% fig2.PaperOrientation = 'landscape';
% fig2.Color = [1 1 1];
% 
% fig2.Position = [0 0 24 12];
% 
% 
% 
% axes1 = axes;
% axes1.Units = 'centimeters';
% axes1.Position = [3 3 8 8];
% axes1.Box = 'on';
% axes1.XColor = 'k';
% axes1.YColor = 'k';
% axes1.Color = 'none';
% axes1.XLabel.String = '1/T (1/K)';
% axes1.YLabel.String = 'ln(k) (mol m^{-2} yr^{-1})';
% hold on
% 
% 
% 
% plOpalA = plot(opalA_temps,opalA_rates,'o');
% plOpalA.MarkerFaceColor = [173, 216, 230]/255; plOpalA.MarkerEdgeColor = plOpalA.MarkerFaceColor/2;
% 
% plBestFit = line([min(opalA_temps),max(opalA_temps)], fitresult.p1*[min(opalA_temps),max(opalA_temps)] + fitresult.p2);
% plBestFit.LineWidth = 1; plBestFit.Color = 'k';
% uistack(plBestFit,'bottom');
% plPredInts = line(xVals,PI);
% plPredInts(1).Color = [0.5 0.5 0.5]; plPredInts(2).Color = [0.5 0.5 0.5];  
% plPredInts(1).LineStyle = '--'; plPredInts(2).LineStyle = '--';
% 
% title('Opal-A dissolution');
% 
% l = legend([plOpalA,plBestFit,plPredInts(1)],{'Experimental data','Best fit line ln(k) = ln(A^*) - E_A/RT','95% prediction intervals'});
% l.Location = 'southwest';
% l.Box = 'off';
% l.Color = 'none';
% 
% 
% % %%%%%%%%%%% quartz:
% 
% axes2 = axes;
% axes2.Units = 'centimeters';
% axes2.Position = [13 3 8 8];
% axes2.Box = 'on';
% axes2.XColor = 'k';
% axes2.YColor = 'k';
% axes2.Color = 'none';
% axes2.XLabel.String = '1/T (1/K)';
% axes2.YLabel.String = 'ln(k) (mol m^{-2} yr^{-1})';
% hold on
% 
% plQtz = plot(quartz_temps,quartz_rates,'o');
% plQtz.MarkerFaceColor = [173, 216, 230]/255; plQtz.MarkerEdgeColor = plQtz.MarkerFaceColor/2;
% 
% plBestFit2 = line([min(quartz_temps),max(quartz_temps)], fitresult2.p1*[min(quartz_temps),max(quartz_temps)] + fitresult2.p2);
% plBestFit2.LineWidth = 1; plBestFit2.Color = 'k';
% uistack(plBestFit2,'bottom');
% plPredInts2 = line(xVals2,PI2);
% plPredInts2(1).Color = [0.5 0.5 0.5]; plPredInts2(2).Color = [0.5 0.5 0.5];  
% plPredInts2(1).LineStyle = '--'; plPredInts2(2).LineStyle = '--';
% 
% title('Quartz dissolution');
% 
% minX = min([axes1.XLim,axes2.XLim]);
% maxX = max([axes1.XLim,axes2.XLim]);
% minY = min([axes1.YLim,axes2.YLim]);
% maxY = max([axes1.YLim,axes2.YLim]);
% 
% axes1.XLim = [minX maxX];
% axes2.XLim = [minX maxX];
% axes1.YLim = [minY maxY];
% axes2.YLim = [minY maxY];
% 
% axes(axes1)
% t1 = text(axes1.XLim(2)*0.98,axes1.YLim(2)*0.95,['E_A = ',num2str(round(Ea(1),1)),char(177) ,num2str(round(errOnEa(1),1)),' kJ mol^{-1}']); t1.HorizontalAlignment = 'right'; t1.VerticalAlignment = 'top';
% t2 = text(axes1.XLim(2)*0.98,axes1.YLim(2)*0.8,['A^* = ',num2str(round(A(1)/10^9,1)),'x10^9 mol m^{-2} yr^{-1}']); t2.HorizontalAlignment = 'right'; t2.VerticalAlignment = 'top';
% t3 = text(axes1.XLim(2)*0.98,axes1.YLim(2)*0.65,['r^2 = ',num2str(round(gof.rsquare,2))]); t3.HorizontalAlignment = 'right'; t3.VerticalAlignment = 'top';
% 
% 
% axes(axes2)
% t4 = text(axes1.XLim(2)*0.98,axes1.YLim(2)*0.95,['E_A = ',num2str(round(Ea(3),1)),char(177) ,num2str(round(errOnEa(3),1)),' kJ mol^{-1}']); t4.HorizontalAlignment = 'right'; t4.VerticalAlignment = 'top';
% t5 = text(axes1.XLim(2)*0.98,axes1.YLim(2)*0.8,['A^* = ',num2str(round(A(3)/10^9,1)),'x10^9 mol m^{-2} yr^{-1}']); t5.HorizontalAlignment = 'right'; t5.VerticalAlignment = 'top';
% t6 = text(axes1.XLim(2)*0.98,axes1.YLim(2)*0.65,['r^2 = ',num2str(round(gof2.rsquare,2))]); t6.HorizontalAlignment = 'right'; t6.VerticalAlignment = 'top';
% 
% 
% axes1.XTickLabel = string( num2cell(axes1.XTick));
% axes2.XTickLabel = string( num2cell(axes2.XTick));

% close(fig2)
cd .., cd ..
