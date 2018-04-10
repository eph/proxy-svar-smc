function vm_plot_irf_bvar(mmodel,prior_type,i_var_instr,p,MP,fflagFEVD,printFig)

graph_opt.font_num = 10;
load(strcat('./results/Result_',char(mmodel),char(i_var_instr),'p_',num2str(p),'_pr_',prior_type,'MP_',num2str(MP),'.mat'))

nv = size(SVAR.LtildeFull,1);
Horizon = size(SVAR.LtildeFull,2);
H = Horizon -1;
nshock = size(SVAR.LtildeFull,4); 
nIV = size(SVAR.i_var_instr,2);
linW = 2;
nshockplot = nIV;
varSelec = SVAR.varSelec;

[nbplt,nr,nc,lr,lc,nstar] = pltorg(length(varSelec));
for jj = 1:nshockplot % Shock
    % Plot IRFs
    figure
    for ii = 1:length(varSelec) % Variable
        subplot(nr,nc,ii)
        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,3,jj)),'k','LineWidth',linW)
        hline(0,'-r')
        hold on
        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,1,jj)),'r--', ...
             'LineWidth',linW)
        plot(0:1:H,squeeze(SVAR.LtildeFull(varSelec(ii),1:Horizon,5,jj)),'r--', ...
             'LineWidth',linW)
        title(SVAR.i_var_str_names(:,varSelec(ii)),'FontSize',graph_opt.font_num,'FontWeight','bold')
        axis([0 H ylim])
        set(gca,'XTick',0:12:H)
        set(gca,'LineWidth',linW)
        grid on
        box off
    end
    if printFig
        print('-dpdf',strcat('./results/IRF_',char(mmodel),char(i_var_instr),'p_',num2str(p),'_priortype_',prior_type, '_MP_',num2str(MP),'.pdf'))

        print('-dpdf',strcat('./results/IRF_',char(mmodel),char(i_var_instr),'p_',num2str(p),'_prtrunc_',prior_type,'_MP_',num2str(MP),'.pdf'))
    end
%     Plot FEVD
    if fflagFEVD ==1
    nvW = length(varSelec);
    [~,nrW,ncW,~,~,~] = pltorg(nvW);
    figure
    for ii = 1:nvW % Variable
        subplot(nrW,ncW,ii)
        plot(0:1:H,squeeze(SVAR.WhFull(varSelec(ii),1:Horizon,3,jj)),'k','LineWidth',linW)
        hold on
        plot(0:1:H,squeeze(SVAR.WhFull(varSelec(ii),1:Horizon,1,jj)),'r--','LineWidth',linW)
        plot(0:1:H,squeeze(SVAR.WhFull(varSelec(ii),1:Horizon,5,jj)),'r--','LineWidth',linW)
        title(SVAR.i_var_str_names(:,varSelec(ii)),'FontSize',graph_opt.font_num,'FontWeight','bold')
        axis([0 H ylim])
        set(gca,'XTick',0:12:H)
        set(gca,'LineWidth',linW)
grid on
box off
    end
    if printFig
        print('-dpdf',strcat('./results/FEVD_',char(mmodel),char(i_var_instr),'p_',num2str(p),'_priortype_',prior_type, '_MP_',num2str(MP),'.pdf'))
    end
    end
end
