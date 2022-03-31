function [X_FL_MPC, T_FL_MPC, U_FL_MPC, X_BAR_FL_MPC, u_FL_MPC, U_FF_FL_MPC] = FL_MPC(p, d, o)
A = [0 1; 0 0];
B = [0; 1];

N = p.MPC.N;
Xf = p.ODE.Xf;
dt = p.MPC.dt;

% MPC size
state_dim = 2;
input_dim = 1;
total_dim = state_dim*N + input_dim*(N-1);

% First FL the system
FL_u = o.FL_u;

[Ad_k,Bd_k,Cd_k] = css2dss(p,dt,A,B,[]);

x = p.ODE.X0;
X_FL_MPC = [];
T_FL_MPC = [];
U_FL_MPC = [];
X_BAR_FL_MPC = [];
U_FF_FL_MPC = [];


yalmip('clear')
vars = sdpvar(total_dim+1,1); % state, inut, and one slack variable
IC = sdpvar(state_dim,1);
Constraints = [];
Constraints = [Constraints vars(state_dim*(N-1)+1:state_dim*N) == Xf];
for i=1:N-1
    Constraints = [Constraints vars(i*state_dim+1:(i+1)*state_dim) == Ad_k*vars((i-1)*state_dim+1:i*state_dim) + Bd_k*vars(state_dim*N+(i-1)*input_dim+1:state_dim*N+i*input_dim)];
end
Constraints = [Constraints kron(diag(ones(1,N)),p.Const.A_in)*vars(1:state_dim*N) <= kron(ones(N,1),p.Const.b_in)+vars(end)];
Constraints = [Constraints vars(state_dim*N+1:total_dim) <= p.Const.u_max];
Constraints = [Constraints vars(state_dim*N+1:total_dim) >= p.Const.u_min];
Constraints = [Constraints vars(1:state_dim) == IC];

Q = eye(total_dim);
Q = [Q zeros(total_dim,1); zeros(1,total_dim) 1e10];
f = zeros(total_dim+1,1);
Objective = 1/2*vars'*Q*vars + f'*vars;

P = optimizer(Constraints,Objective,sdpsettings('solver','mosek','verbose',0),...
    IC,vars);

for iter = 1:p.ODE.tspan(end)/dt
    disp(iter)
    
    x0 = x(end,:)';
    %     Constraints(end) = vars(1:state_dim) == x0';
    
    [sol, diagnostics] = P(x0);
    if diagnostics ~= 0
        error('Issue with Mosek in FL+MPC');
    end
    
    MPC_sol = sol;
    
    t_FL_MPC = 0:dt:dt*N;
    x_FL_MPC = [MPC_sol(1:state_dim:state_dim*N) MPC_sol(2:state_dim:state_dim*N)];
    u_lin_k = MPC_sol(state_dim*N+1:total_dim);
    slack = MPC_sol(end);
    
    [t,x] = ode45(@(t,x) d.f_func_w(x(1),x(2),t) + d.g_func_w(x(1),x(2),t)*(FL_u(x(1),x(2),u_lin_k(1))),[0 dt],x0);
    
    clear u_FL_MPC;
    for i=1:length(x)
        u_FL_MPC(i) = FL_u(x(i,1),x(i,2),u_lin_k(1));
    end
    
    X_BAR_FL_MPC = [X_BAR_FL_MPC; x_FL_MPC(2,:)];
    
    if isempty(T_FL_MPC)
        X_FL_MPC = [X_FL_MPC; x];
        T_FL_MPC = [T_FL_MPC; t];
        U_FL_MPC = [U_FL_MPC; u_FL_MPC'];
        U_FF_FL_MPC = [U_FF_FL_MPC; u_lin_k(1)*ones(size(t))];
    else
        T_FL_MPC = [T_FL_MPC; t(2:end)+T_FL_MPC(end)];
        X_FL_MPC = [X_FL_MPC; x(2:end,:)];
        U_FL_MPC = [U_FL_MPC; u_FL_MPC(2:end)'];
        U_FF_FL_MPC = [U_FF_FL_MPC; u_lin_k(1)*ones(size(t(2:end)))];
    end
end
end