function updateProgressPlots(z_real,massStruct,ratioStruct,solubilities,dt,porespace,a1a, a1b, a1c, a2, a3,a4, i, sedAmt)
        %updateProgressPlots(z_real,massStruct,ratioStruct,solubilities,dt,porespace,a1a, a1b, a1c, a2, a3, ii)
 
% unpack structure...

%z_real =            varsForPlotting.z_real;
porewater_Si_conc = massStruct.porewater_Si ./ porespace;
SolQtz =            solubilities.quartz;
SolCT=              solubilities.opalCT;
SolASi =            solubilities.opalA;
opalA_Si =          massStruct.opalA_Si;
qtz_Si =            massStruct.quartz_Si;
opalCT_Si =         massStruct.opalCT_Si;
d18O_porewater =    ratioStruct.porewater_d18O;
d18O_opal =         ratioStruct.opalA_d18O;
d18O_CT =           ratioStruct.opalCT_d18O;
d18O_qtz =          ratioStruct.quartz_d18O;

d30Si_porewater =    ratioStruct.porewater_d30Si;
d30Si_opal =         ratioStruct.opalA_d30Si;
d30Si_CT =           ratioStruct.opalCT_d30Si;
d30Si_qtz =          ratioStruct.quartz_d30Si;

%d18O_qtz_inst =     varsForPlotting.d18O_qtz_inst;
%varsForPlotting.mean_spacing


            % Plot1, porewater Si, + solubilities and temperature on secondary axes for info.
            axes(a1a); 
            plot(z_real(1:i),porewater_Si_conc(1:i)*10^3,'k:','LineWidth',2); a1a.Color = 'none';
            hold on
            plot(z_real(1:i),SolQtz(1:i)*10^6,'r--'); text(z_real(end-10),SolQtz(end)*10^6,'Quartz solubility','Color','r','HorizontalAlignment','right');
            plot(z_real(1:i),SolCT(1:i)*10^6,'b--'); text(z_real(end-10),SolCT(end)*10^6,'OpalCT solubility','Color','b','HorizontalAlignment','right');
            plot(z_real(1:i),SolASi(1:i)*10^6,'k-'); text(z_real(end-10),SolASi(end)*10^6,'OpalA solubility','Color','k','HorizontalAlignment','right');
            hold off,
            a1a.XLim = a1b.XLim;
            
            % pl = plot(a1c,IODP_table.Depth_mbsf_,IODP_table.Silica_H4SiO4__uM_/1000,'.'); pl.MarkerEdgeColor = [0.7 0.7 0.7];
            a1c.YLim = a1a.YLim; a1c.XLim = a1b.XLim;
            a1c.Color = 'none';  a1c.XTick = [];  a1c.YTick = [];
            
            a1a.YLabel.String = 'Porewater Si conc (umolar)'; a1a.XLabel.String = 'Depth (m)'; drawnow;
            
            % Plot2
            axes(a2);
            plot(z_real(1:i),opalA_Si(1:i)*60/sedAmt*100,'b-o');  hold on
            %scatter(z_real(1:i),immatureCT(1:i),5,mean_spacing(1:i),'filled')
            plot(z_real(1:i),qtz_Si(1:i)*60/sedAmt*100,'k:');  title([num2str(i*dt/1000),' Myr'])
            plot(z_real(1:i),opalCT_Si(1:i)*60/sedAmt*100,'g-'), hold off
            legend({'opalA','quartz','opalCT'})
            a2.YLabel.String = 'sediment concentration (wt% dry sed)'; a2.XLabel.String = 'Depth (m)';  
            a2.XLim = a1b.XLim;
            a2.YLim(1) = 0; drawnow;
            
            % Plot3
            axes(a3);
            plot(z_real(1:i), d18O_porewater(1:i),  'k-');      hold on
            plot(z_real(1:i), d18O_opal(1:i),       'b--');
            plot(z_real(1:i), d18O_CT(1:i),         'r.');
            plot(z_real(1:i), d18O_qtz(1:i),        'g-');
            %plot(z_real(1:i-1), d18O_qtz_inst(1:i-1),   'k.')
            hold off
            a3.YLim = [-5 45];
            legend({'porewater','opalA','opalCT','cumulative quartz'});%,'instantaneous quartz'})
            a3.YLabel.String = 'Porewater d18O (permil)'; a3.XLabel.String = 'Depth (m)'; 
            a3.XLim = a1b.XLim;
            drawnow;

            % Plot3
            axes(a4);
            plot(z_real(1:i), d30Si_porewater(1:i),  'k-');      hold on
            plot(z_real(1:i), d30Si_opal(1:i),       'b--');
            plot(z_real(1:i), d30Si_CT(1:i),         'r.');
            plot(z_real(1:i), d30Si_qtz(1:i),        'g-');
            %plot(z_real(1:i-1), d30Si_qtz_inst(1:i-1),   'k.')
            hold off
            a4.YLim = [-5 5];
            legend({'porewater','opalA','opalCT','cumulative quartz'});%,'instantaneous quartz'})
            a4.YLabel.String = 'Porewater d30Si (permil)'; a4.XLabel.String = 'Depth (m)'; 
            a4.XLim = a1b.XLim;
            drawnow;

end