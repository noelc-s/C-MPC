
figure(1);
set(gcf,'Position',[1000 0 800 800]);
clf;
% subplot(1,2,1)
hold on;
col = [136 34 85;
       136 204 238;
       22 119 51;
       51 34 136]/255; % color blind friendly palette
plot(X_Lin_MPC_CLF(:,1),X_Lin_MPC_CLF(:,2),'linewidth',2,'color',col(4,:))
plot(XD_Lin_MPC(:,1),XD_Lin_MPC(:,2),'--','linewidth',2,'color',col(4,:))
scatter(X_K_MPC_CLF(:,1),X_K_MPC_CLF(:,2),50,col(4,:),'filled','linewidth',2)
yline(p.Const.b_in(3),'linewidth',2)
yline(-p.Const.b_in(4),'linewidth',2)
xline(p.Const.b_in(1),'linewidth',2)
xline(-p.Const.b_in(2),'linewidth',2)
axis([-3 0.5 -2.1 .8])
% axis([-1.5 -0.4 0.3 .6])
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'linewidth',2)
% subplot(1,2,2)
% hold on;
% plot(X_Lin_MPC_CLF(:,1),X_Lin_MPC_CLF(:,2),'linewidth',2,'color',col(4,:))
% scatter(X_K_MPC_CLF(:,1),X_K_MPC_CLF(:,2),50,col(4,:),'filled','linewidth',2)
% yline(p.Const.b_in(3),'linewidth',2)
% yline(-p.Const.b_in(4),'linewidth',2)
% xline(p.Const.b_in(1),'linewidth',2)
% xline(-p.Const.b_in(2),'linewidth',2)
% % axis([-3 0.5 -2.1 .8])
% axis([-1.5 -0.4 0.3 .6])
% xlabel('$x_1$','interpreter','latex')
% ylabel('$x_2$','interpreter','latex')
% set(gca,'TickLabelInterpreter', 'latex');
% set(gca,'FontSize',20)
% set(gca,'linewidth',2)
