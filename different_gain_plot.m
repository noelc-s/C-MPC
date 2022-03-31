figure(1);
set(gcf,'Position',[1000 0 1000 500]);
clf;
col = parula(16);
for k = 1:2:15
    subplot(1,2,1)
    hold on;
    plot(X_Lin_MPC_CLF{k}(:,1),X_Lin_MPC_CLF{k}(:,2),'linewidth',2,'color',col(k,:))
    scatter(X_K_MPC_CLF{k}(:,1),X_K_MPC_CLF{k}(:,2),50,col(k,:),'filled','linewidth',2)
    yline(p.Const.b_in(3),'linewidth',2)
    yline(-p.Const.b_in(4),'linewidth',2)
    xline(p.Const.b_in(1),'linewidth',2)
    xline(-p.Const.b_in(2),'linewidth',2)
%     axis([0 1.3 -1 1.1])
    axis([0.9 1.3 0 1.1])
    xlabel('$\theta$','interpreter','latex')
    ylabel('$\dot\theta$','interpreter','latex')
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',20)
    set(gca,'linewidth',2)
    subplot(1,2,2)
    hold on;
    plot(T_Lin_MPC_CLF{k},U_Lin_MPC_CLF{k},'linewidth',2,'color',col(k,:))
    plot(T_Lin_MPC_CLF{k},U_FF_MPC_CLF{k},'--','linewidth',2,'color',col(k,:))
    yline(p.Const.u_max,'linewidth',2)
    yline(p.Const.u_min,'linewidth',2)
    axis([0 1.5 p.Const.u_min-0.5 1])
    xlabel('t','interpreter','latex')
    ylabel('u','interpreter','latex')
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',20)
    set(gca,'linewidth',2)
end