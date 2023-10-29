%% Initialize

init();

p = params();
[d, o] = dynamics_and_outputs(p);

%% FL
v = -p.FL.alpha_FL*[o.y; o.Lfy]; % auxiliary controller
u = o.LgLfy\(-o.Lf2y + v); % nonlinear controller
FL = o.FL_u;

% Solve ODE
[t,x] = ode45(@(t,x) d.f_func_w(x(1),x(2),t) + d.g_func_w(x(1),x(2),t)*FL(x(1),x(2),-p.FL.alpha_FL*x),p.ODE.tspan,p.ODE.X0);

% Gather ODE output
t_FL = t;
x_FL = x;
u_FL = FL(x(:,1),x(:,2),x*p.FL.alpha_FL');

%% CLF

V = @(x) o.eta(x(1),x(2))'*p.CLF.P_lyap*o.eta(x(1),x(2));
LFV = @(x) o.eta(x(1),x(2))'*(p.CLF.F'*p.CLF.P_lyap + p.CLF.P_lyap*p.CLF.F)*o.eta(x(1),x(2));
LGV = @(x) 2*o.eta(x(1),x(2))'*p.CLF.P_lyap*p.CLF.G;

options = optimset('Display', 'off');
% FL-CLF-QP_lyap
v = @(x) quadprog(1,0, LGV(x), -p.CLF.gamma*V(x)-LFV(x),[],[],[],[],[],options);
u = @(x) o.LgLfy_func(x(1),x(2))\(-o.Lf2y_func(x(1),x(2)) + v(x));
[t,x] = ode45(@(t,x) d.f_func_w(x(1),x(2),t) + d.g_func_w(x(1),x(2),t)*u(x),p.ODE.tspan,p.ODE.X0);
t_CLF = t;
x_CLF = x;
u_CLF = zeros(size(t));
slack = zeros(size(t));
for i=1:length(t)
    u_CLF(i) = u(x(i,:));
end

%% FL, then MPC on linear system
[X_FL_MPC, T_FL_MPC, U_FL_MPC, X_BAR_FL_MPC, u_FL_MPC, U_FF_FL_MPC] = FL_MPC(p, d, o);

%% MPC-CLF proposed approach
% p.low_level = 'None';
p.low_level = 'CLF';
[X_Lin_MPC_CLF, T_Lin_MPC_CLF, U_Lin_MPC_CLF, X_K_MPC_CLF, XD_Lin_MPC, u_Lin_MPC_CLF, U_FF_MPC_CLF] = MPC_Bez(p, d, o, T_FL_MPC, X_FL_MPC, U_FL_MPC);

%% Plotting
proposed_plot();