function a = solubilitySuppression(varargin)

% PJF September 2021. Takes literature datasets from Dixit et al,
% Gallinari et al, McManus et al, and van Cappellen and Qui, that have
% published detrital mineral abundance and opal (mostly/entirely biogenic
% silica) solubility. Uses this data to define a relationship between the
% 'suppression' of opal solubility and the amount of detrital
% aluminosilicates present. Essentially it is just fitting a quadratic
% through all the data, where each individual dataset is normalised to the
% maximum/mineral-free solubility (to account for differences in e.g.
% intrinsic opal reactivity, surface area, temperature of experiments,
% etc., that might all influence solubility). Choice of a quadratic
% function is largely arbitrary, but it fits the data well.
%

if ~isempty(varargin)
    if varargin{1} == 1
        plotting = 1;
    else
        plotting = 0;
    end
else
    plotting = 0;
end

% The mixtures represented by the open symbols were prepared by keeping the total mass of solid constant at 100 mg

cd auxFiles/inputData, X = importdata('DixitKaolinite.txt');
KaolOpalRatio = X.data(:,1);        % Detrital minerals percent / opal percent. From Dixit et al. 2001
KaolSolubility = X.data(:,2);       % solubility in batch reactor at 20 deg C
KaolOpalRatio(end+1) = 0;
KaolSolubility(end+1) = 1020;
OpalPC = 1 ./ (KaolOpalRatio + 1);
KaolPC = 1 - OpalPC;
KaolRelSol = KaolSolubility/1020.8 * -1 + 1;


X = importdata('DixitBasalt.txt');  % From Dixit et al. 2001
BasaltOpalRatio = X.data(:,1);
BasaltSolubility = X.data(:,2);
BasaltOpalRatio(end+1) = 0;
BasaltSolubility(end+1) = 1020;
OpalPC_bas = 1 ./ (BasaltOpalRatio + 1);
BasaltPC = 1 - OpalPC_bas;
BasaltRelSol = BasaltSolubility/1020 * -1 + 1;


X = importdata('Thalassiosira.txt'); % This data from Gallinari et al.
KaolThalassRatio = X.data(:,1);
ThalassSolubility = X.data(:,2);
OpalPC_Thalass = 1 ./ (KaolThalassRatio + 1);
KaolPC_Thalass = 1 - OpalPC_Thalass;
ThalassRelSol = ThalassSolubility/ThalassSolubility(1) * -1 + 1;


X = importdata('Cylindrotheca.txt'); % This data from Gallinari et al.
KaolCylindroRatio = X.data(:,1);
CylindroSolubility = X.data(:,2);
OpalPC_Cylindro = 1 ./ (KaolCylindroRatio + 1);
KaolPC_Cylindro = 1 - OpalPC_Cylindro;
CylindroRelSol = CylindroSolubility/CylindroSolubility(1) * -1 + 1;

X = importdata('VC and Q.txt');
VCQ_Ratio = X.data(:,1);
VCQ_Solubility = X.data(:,2);
VCQ_opalPC = 1 ./ (VCQ_Ratio + 1);
VCQ_detritalPC = 1 - VCQ_opalPC;
VCQ_RelSol = VCQ_Solubility / VCQ_Solubility(1) * -1 + 1;


X = importdata('McManus from Gallinari.txt');
McG_Ratio = X.data(:,1);
McG_Solubility = X.data(:,2);
McG_opalPC = 1./ (McG_Ratio + 1);
McG_detPC = 1 - McG_opalPC;
McG_RelSol = McG_Solubility / 600 * -1 + 1; % Eyeballed for now from Gallinari graph - need to have the raw data from MaManus, but no det mins in there. Gallinari cite a Rageauneau et al report...?

allDataX = [KaolPC; BasaltPC; KaolPC_Thalass; KaolPC_Cylindro ;VCQ_detritalPC ;McG_detPC ];
allDataY = [KaolRelSol; BasaltRelSol; ThalassRelSol; CylindroRelSol ; VCQ_RelSol; McG_RelSol];

[xData, yData] = prepareCurveData( allDataX, allDataY );
ft = fittype( 'a*x^2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0;
[fitresult, gof] = fit( xData, yData, ft, opts );

a = fitresult.a;

cd .., cd ..

if plotting == 1
    
    
    [fig,ax] = pjfSquarePlot;
    fig.Name = 'Degree of solubility suppression caused by detrital minerals';
    
    ax.XLim = [0 1];
    ax.YLim = [0 1.05];
    
    box on
    axis square
    hold on
    
    p1 = plot(KaolPC,KaolRelSol,'ro');
    p2 = plot(BasaltPC,BasaltRelSol,'kd');
    p3 = plot(KaolPC_Thalass,ThalassRelSol,'gs');
    p4 = plot(KaolPC_Cylindro,CylindroRelSol,'bp');
    
    p5 = plot(VCQ_detritalPC,VCQ_RelSol,'bd');
    p6 = plot(McG_detPC,McG_RelSol,'rs');
    
    p1.MarkerFaceColor = [1 0 0]; p1.MarkerEdgeColor = p1.MarkerFaceColor * 0.6; p1.Marker = 'o';
    p2.MarkerFaceColor = [0 1 0]; p2.MarkerEdgeColor = p2.MarkerFaceColor * 0.6; p2.Marker = 's';
    p3.MarkerFaceColor = [0 0 1]; p3.MarkerEdgeColor = p3.MarkerFaceColor * 0.6; p3.Marker = 'd';
    p4.MarkerFaceColor = [0.4 0.4 0.4]; p4.MarkerEdgeColor = p4.MarkerFaceColor * 0.6; p4.Marker = 'p';
    p5.MarkerFaceColor = [1 1 0]; p5.MarkerEdgeColor = p5.MarkerFaceColor * 0.6; p5.Marker = '^';
    p6.MarkerFaceColor = [0 1 1]; p6.MarkerEdgeColor = p6.MarkerFaceColor * 0.6; p6.Marker = '>';
    
    ax.XLabel.String = 'Fraction detrital minerals';
    ax.YLabel.String = 'OpalA solubility suppression';
    
    pFITLine1 = plot(0:0.01:1,fitresult.a*(0:0.01:1).^2);
    
    l = legend([p1 p2 p3 p4 p5 p6, pFITLine1],{'Biosiliceous Ooze + Kaolinite^1',...
        'Biosiliceous Ooze + Basalt^1',...
        'T. pseudonanna + Kaolinite^2',...
        'C. fusiformis + Kaolinite^2',...
        'Southern Ocean porewater asymptotes^3',...
        'Equatorial Pacific porewater asymptotes^{2,4}',...
        ['Fit to all data, r^2 = ',num2str(round(gof.rsquare,2))]});
    l.Location = 'northwest';
    l.Box = 'off';
    
    
    
end

