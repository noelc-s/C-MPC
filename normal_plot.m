
figure(1);
set(gcf,'Position',[1000 0 1000 400]);
clf;
subplot(1,2,1)
hold on;
col = [136 34 85;
       136 204 238;
       22 119 51;
       51 34 136]/255; % color blind friendly palette
plot(x_CLF(:,1),x_CLF(:,2),'linewidth',2,'color',col(1,:))
% plot(X_FL_MPC(:,1), X_FL_MPC(:,2),'linewidth',2,'color',col(2,:))
% scatter(X_BAR_FL_MPC(:,1),X_BAR_FL_MPC(:,2),50,col(2,:),'filled','linewidth',2)
plot(X_MPC_FL(:,1), X_MPC_FL(:,2),'linewidth',2,'color',col(3,:))
scatter(X_BAR_MPC_FL(:,1),X_BAR_MPC_FL(:,2),50,col(3,:),'filled','linewidth',2)
plot(X_Lin_MPC_CLF(:,1),X_Lin_MPC_CLF(:,2),'linewidth',2,'color',col(4,:))
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
% legend({'CLF','MPC On FL','','MPC on Linearization with No Low Level','',' MPC on Linearization, CLF Low Level','','','',''},'interpreter','latex')
% legend({'CLF','MPC on Linearization with No Low Level','',' MPC on Linearization, CLF Low Level','','','',''},'interpreter','latex')
subplot(1,2,2)
hold on;
plot(x_CLF(:,1),x_CLF(:,2),'linewidth',2,'color',col(1,:))
% plot(X_FL_MPC(:,1), X_FL_MPC(:,2),'linewidth',2,'color',col(2,:))
% scatter(X_BAR_FL_MPC(:,1),X_BAR_FL_MPC(:,2),50,col(2,:),'filled','linewidth',2)
plot(X_MPC_FL(:,1), X_MPC_FL(:,2),'linewidth',2,'color',col(3,:))
scatter(X_BAR_MPC_FL(:,1),X_BAR_MPC_FL(:,2),50,col(3,:),'filled','linewidth',2)
plot(X_Lin_MPC_CLF(:,1),X_Lin_MPC_CLF(:,2),'linewidth',2,'color',col(4,:))
scatter(X_K_MPC_CLF(:,1),X_K_MPC_CLF(:,2),50,col(4,:),'filled','linewidth',2)
yline(p.Const.b_in(3),'linewidth',2)
yline(-p.Const.b_in(4),'linewidth',2)
xline(p.Const.b_in(1),'linewidth',2)
xline(-p.Const.b_in(2),'linewidth',2)
% axis([-3 0.5 -2.1 .8])
axis([-1.5 -0.4 0.3 .6])
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'linewidth',2)





% hold on;
% plot(t_CLF,u_CLF,'linewidth',2,'color',col(1,:))
% % plot(T_FL_MPC,U_FL_MPC,'linewidth',2,'color',col(2,:))
% % plot(T_FL_MPC,U_FF_FL_MPC,'--','linewidth',2,'color',col(2,:))
% plot(T_MPC_FL,U_MPC_FL,'linewidth',2,'color',col(3,:))
% plot(T_Lin_MPC_CLF,U_Lin_MPC_CLF,'linewidth',2,'color',col(4,:))
% plot(T_Lin_MPC_CLF,U_FF_MPC_CLF,'--','linewidth',2,'color',col(4,:))
% yline(p.Const.u_max,'linewidth',2)
% yline(p.Const.u_min,'linewidth',2)
% axis([0 floor(p.ODE.tspan(end)*2/3) -1 p.Const.u_max+0.2])
% % axis([0 0.5 1 18])
% xlabel('t','interpreter','latex')
% ylabel('u','interpreter','latex')
% set(gca,'TickLabelInterpreter', 'latex');
% set(gca,'FontSize',20)
% set(gca,'linewidth',2)