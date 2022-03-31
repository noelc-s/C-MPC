%% Initialize

init();

p = params_slow();
% p = params();
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

% FL-CLF-QP_lyap with input constraints
% calF = @(x) o.Lf2y_func(x(1),x(2));
% calG = @(x) o.LgLfy_func(x(1),x(2));
% v = @(x) quadprog([calG(x)'*calG(x) 0; 0 1e10],[calF(x)'*calG(x); 0], [LGV(x)*calG(x) -1; 1 0; -1 0], [-p.CLF.gamma*V(x)-LFV(x)-LGV(x)*calF(x); p.Const.u_max; -p.Const.u_min],[],[],[],[],[],options);
% u = @(x) subsref(v(x),struct('type','()','subs',{{1}}));
% delta = @(x) subsref(v(x),struct('type','()','subs',{{2}}));

[t,x] = ode45(@(t,x) d.f_func_w(x(1),x(2),t) + d.g_func_w(x(1),x(2),t)*u(x),p.ODE.tspan,p.ODE.X0);
t_CLF = t;
x_CLF = x;
u_CLF = zeros(size(t));
slack = zeros(size(t));
for i=1:length(t)
    u_CLF(i) = u(x(i,:));
%     slack(i) = delta(x(i,:));
end

%% FL, then MPC on linear system
[X_FL_MPC, T_FL_MPC, U_FL_MPC, X_BAR_FL_MPC, u_FL_MPC, U_FF_FL_MPC] = FL_MPC(p, d, o);

%% MPC on linearized system
disp('MPC on linearized system');
p.low_level = 'None';
[X_MPC_FL, T_MPC_FL, U_MPC_FL, X_BAR_MPC_FL, u_MPC_FL, U_FF_MPC_FL] = MPC_FL(p, d, o,T_FL_MPC, X_FL_MPC, U_FL_MPC);

%% MPC-CLF proposed approach
% p.low_level = 'None';
p.low_level = 'CLF';
[X_Lin_MPC_CLF, T_Lin_MPC_CLF, U_Lin_MPC_CLF, X_K_MPC_CLF, u_Lin_MPC_CLF, U_FF_MPC_CLF] = MPC_Bez(p, d, o, T_FL_MPC, X_FL_MPC, U_FL_MPC);

%% Plotting
slow_plot();
% normal_plot();