function [X_Lin_MPC, T_Lin_MPC, U_Lin_MPC, X_K_MPC_CLF, u_Lin_MPC, U_FF_MPC_CLF] = MPC_FL(p, dyn, o, T_FL_MPC, X_FL_MPC, U_FL_MPC)
disp('Setting up standard approach');

% TODO: Q1) need to use a feasible (x,u) to start, which the following is not.
t_bar = T_FL_MPC;
x_bar = X_FL_MPC;
u_bar = U_FL_MPC;

x = p.ODE.X0;
X_Lin_MPC = [];
T_Lin_MPC = [];
U_Lin_MPC = [];
X_K_MPC_CLF = [];
U_FF_MPC_CLF = [];

N = p.MPC.N;
Xf = p.ODE.Xf;
dt = p.MPC.dt;
state_stage_cost = p.MPC_CLF.state_stage_cost;
input_stage_cost = p.MPC_CLF.input_stage_cost;
aux_stage_cost = p.MPC_CLF.aux_stage_cost;
A_in = p.Const.A_in;
b_in = p.Const.b_in;
P_lyap = p.CLF.P_lyap;

% MPC size
state_dim = 2;
state_size_without_terminal = state_dim*(N-1);
state_size_with_terminal = state_dim*N;
state_end_index=state_dim*N;
input_dim = 1;
input_size = input_dim*(N-1);
input_end_index=state_end_index+input_dim*(N-1);
auxiliary_dim = 2;
auxiliary_size = auxiliary_dim*(N-1);
auxiliary_end_index=input_end_index+auxiliary_dim*(N-1);
total_dim = state_dim*N + input_dim*(N-1) + auxiliary_dim*(N-1);

% Clf redefinition
eta_point = o.eta;
V_point = @(x) eta_point(x(1),x(2))'*P_lyap*eta_point(x(1),x(2));

yalmip('clear')
% P_lyapre-size
A = cell(N-1);
B = cell(N-1);
C = cell(N-1);
Ad_k = zeros(state_dim, state_size_without_terminal);
Bd_k = zeros(state_dim, N-1);
Cd_k = zeros(state_dim, N-1);
x_lin_k = zeros(N-1,2);
u_lin_k = zeros(N-1,1);
Ad_km1 = zeros(state_dim, state_size_without_terminal);
Bd_km1 = zeros(state_dim, N-1);
Cd_km1 = zeros(state_dim, N-1);
x_lin_km1 = zeros(N-1,2);
u_lin_km1 = zeros(N-1,1);

% Optimizer variables
states = sdpvar(2,N);
inputs = sdpvar(input_size,1);
slack = sdpvar(1);
Ad_ = sdpvar(state_dim, state_size_without_terminal);
Bd_ = sdpvar(state_dim, N-1);
Cd_ = sdpvar(state_dim, N-1);
x_lin_ = sdpvar(N-1,2);
u_lin_ = sdpvar(N-1,1);
x0_ = sdpvar(state_dim,1);

% adding a slack variable for state constraints
Constraints = [];

% clear u_k s
u_k = sdpvar(N-1,1);

for i = 1:N-1
    % Decision variables
    Constraints = [Constraints u_k(i,1) == inputs((i-1)*input_dim+1:i*input_dim)];
    
    % Dynamics
    Constraints = [Constraints states(:,i+1) == Ad_(1:state_dim,(i-1)*state_dim+1:(i)*state_dim)*states(:,i)...
        + Bd_(1:state_dim, i)*u_k(i,1) + Cd_(1:state_dim, i)];
    %     Constraints = [Constraints norm(sqrtm(M_k)*s(:,i) + 1/2*inv(sqrtm(M_k))*N_k_(:,i),2) <= (1/4*N_k_(:,i)'*inv(M_k)*N_k_(:,i)-(u_max- Gamma_k_(i)))^(1/2)];
    %     Constraints = [Constraints s(:,i)'*M_k*s(:,i) + N_k_(:,i)'*s(:,i) + Gamma_k_(i) <= u_max];
    
    % b_k in X \ominus E
    Constraints = [Constraints u_k(i,1) <= p.Const.u_max];%+slack];
    Constraints = [Constraints u_k(i,1) >= p.Const.u_min];%+slack];
    Constraints = [Constraints A_in*states(:,i) <= b_in+slack];
    %     Constraints = [Constraints norm(s(:,i),2) <= 1];
end

% Terminal Set
Constraints = [Constraints states(:,end) == Xf];

% Add this for ISS tube, which is not considered yet.
% Constraints = [Constraints V_point(states(1:state_dim)-x0_) <= gamma_MPCFL];
Constraints = [Constraints states(:,1) == x0_];

Q = state_stage_cost*eye(state_dim*N);
R = input_stage_cost*eye(input_dim*(N-1));
W = 1e9;
% Q = blkdiag(Q,R,Aux,1e9);
% f = zeros(total_dim+1,1);
Objective = 1/2*(states(:)'*Q*states(:) + inputs'*R*inputs + slack'*W*slack);

% SP_lyapQ_epsilon = 0.1;
% SQP_lyapconverged = false;

P = optimizer(Constraints,Objective,sdpsettings('solver','mosek','verbose',0),...
    {Ad_,Bd_,Cd_,x_lin_,u_lin_,x0_},...
    {states,inputs,slack,u_k});

%%% Get linearizations at origin for backup plan
x_F(i,:) = [0 0];
u_F(i,:) = 0;
i = 1;
A_origin = dyn.Df_func_model(x_F(i,1), x_F(i,2))+dyn.Dg_func_model(x_F(i,1), x_F(i,2))*u_F(i,:);
B_origin = dyn.g_func_model(x_F(i,1),x_F(i,2));
C_origin = dyn.f_func_model(x_F(i,1),x_F(i,2)) + dyn.g_func_model(x_F(i,1),x_F(i,2))*u_F(i,:) - (A_origin*x_F(i,:)' + B_origin*u_F(i,:));
[Ad_,Bd_,Cd_] = css2dss(p,dt,A_origin,B_origin,C_origin);
Ad_at_origin = Ad_;
Bd_at_origin = Bd_;
Cd_at_origin = Cd_;

for iter = 1:p.ODE.tspan(end)/dt
    disp(iter)
    x0 = x(end,:)';
    %     while(~SQP_lyapconverged)
    for i = 1:N-1
        x_lin_k(i,:) = interp1(t_bar, x_bar, dt*(i-1));
        u_lin_k(i,:) = interp1(t_bar, u_bar, dt*(i-1));
        A{i} = dyn.Df_func_model(x_lin_k(i,1), x_lin_k(i,2))+dyn.Dg_func_model(x_lin_k(i,1), x_lin_k(i,2))*u_lin_k(i,:);
        B{i} = dyn.g_func_model(x_lin_k(i,1),x_lin_k(i,2));
        C{i} = dyn.f_func_model(x_lin_k(i,1),x_lin_k(i,2)) + dyn.g_func_model(x_lin_k(i,1),x_lin_k(i,2))*u_lin_k(i,:) - (A{i}*x_lin_k(i,:)' + B{i}*u_lin_k(i,:));
        [Ad_,Bd_,Cd_] = css2dss(p,dt,A{i},B{i},C{i});
        Ad_k(1:state_dim,(i-1)*state_dim+1:(i)*state_dim) = Ad_;
        Bd_k(1:state_dim, i) = Bd_;
        Cd_k(1:state_dim, i) = Cd_;
    end
    
    % Constrain acceleration of starting point
    %         inputs = vars(state_end_index+1:input_end_index);
    %         clear f_ g_ f_eval g_eval;
    %         for i = 1:N-1
    %             pos1 = vars((i-1)*state_dim+1);
    %             vel1 = vars((i-1)*state_dim+2);
    %             pos4 = vars((i)*state_dim+1);
    %             vel4 = vars((i)*state_dim+2);
    %             % Defined:
    %             pos2 = dt/3*vel1+pos1;
    %             pos3 = -dt/3*vel4+pos4;
    %             vel2 = H(2,:)*[pos1 pos2 pos3 pos4]';
    %             vel3 = H(3,:)*[pos1 pos2 pos3 pos4]';
    %             p_dd_0(i) = H(1,:)*[vel1 vel2 vel3 vel4]';
    %             f_eval = f_func_model(x_lin(i,1),x_lin(i,2));
    %             f_(i) = f_eval(2);
    %             g_eval = g_func_model(x_lin(i,1),x_lin(i,2));
    %             g_(i) = g_eval(2);
    %         end
    
    %         diagnostics = optimize(Constraints, 1/2*vars'*Q*vars + f'*vars + ...
    %             0*norm(p_dd_0' - f_'-g_'.*inputs,2).^2,opts);
    [sol, diagnostics,d1,d2,d3,d4] = P({Ad_k,Bd_k,Cd_k,x_lin_k,u_lin_k,x0});
    if iter == 1 & diagnostics ~= 0
        error('Issue with Mosek in proposed');
    elseif diagnostics ~= 0
        warning('Falling back on previous linearization')
        Ad_km1 = [Ad_km1(:,3:end) Ad_at_origin];
        Bd_km1 = [Bd_km1(:,2:end) Bd_at_origin];
        Cd_km1 = [Cd_km1(:,2:end) Cd_at_origin];
        % TODO: assumes terminal point is unforced equilibrium
        x_lin_km1 = [x_lin_km1(2:end,:); 0 0];
        u_lin_km1 = [u_lin_km1(2:end,:); 0];
        [sol, diagnostics,d1,d2,d3,d4] = P({Ad_km1,Bd_km1,Cd_km1,x_lin_km1,u_lin_km1,x0});
    else
        Ad_km1 = Ad_k;
        Bd_km1 = Bd_k;
        Cd_km1 = Cd_k;
        x_lin_km1 = x_lin_k;
        u_lin_km1 = u_lin_k;
    end
    
    %     MPC_sol = value(vars);
    
    t_FL_MPC = 0:dt:dt*(N-1);
    %     x_FL_MPC = [MPC_sol(1:2:state_end_index) MPC_sol(2:2:state_end_index)];
    %     u_lin = MPC_sol(state_end_index+1:input_end_index);
    %     aux_vars = MPC_sol(input_end_index+1:auxiliary_end_index);
    %     slack = MPC_sol(end);
    x_FL_MPC = sol{1};
    u_lin_k = sol{2};
    aux_vars = sol{3};
    slack = sol{4};
    
    %         if norm(x_FL_MPC(1:end-1,:)-x_lin,2) < SP_lyapQ_epsilon
    %             SQP_lyapconverged = true;
    %         else
    %             t_bar = t_FL_MPC;
    %             x_bar = x_FL_MPC;
    %             u_bar = [u_lin; 0];
    %         end
    %     end
    %     SQP_lyapconverged = false;
    
    %%%%%%%%%%%%%%%%%%%%% Which model to use to reconstruct the continuous
    %%%%%%%%%%%%%%%%%%%%% time trajectory to track? %%%%%%%%%%%%%%%%%%%%%%
    %%% Bezier
    
    [t,x] = ode45(@(t,x) dyn.f_func_w(x(1),x(2),t) + dyn.g_func_w(x(1),x(2),t)*(u_lin_k(1)),[0 dt],x0); % no low level
    clear u_Lin_MPC;
    for i=1:length(x)
        switch p.low_level
            case 'None'
                u_Lin_MPC(i) = u_lin_k(1);
            case 'CLF'
                u_Lin_MPC(i) = FL_CLF1(t(i),x(i,:));
                %         u_Lin_MPC(i) = FL_CLF2(t(i),x(i,:));
        end
    end
    
    if isempty(T_Lin_MPC)
        T_Lin_MPC = [T_Lin_MPC; t];
        X_Lin_MPC = [X_Lin_MPC; x];
        U_Lin_MPC = [U_Lin_MPC; u_Lin_MPC'];
        U_FF_MPC_CLF = [U_FF_MPC_CLF; u_lin_k(1)*ones(size(t))];
    else
        T_Lin_MPC = [T_Lin_MPC; t(2:end)+T_Lin_MPC(end)];
        X_Lin_MPC = [X_Lin_MPC; x(2:end,:)];
        U_Lin_MPC = [U_Lin_MPC; u_Lin_MPC(2:end)'];
        U_FF_MPC_CLF = [U_FF_MPC_CLF; u_lin_k(1)*ones(size(t(2:end)))];
    end
    
    X_K_MPC_CLF = [X_K_MPC_CLF; x_FL_MPC(:,2)'];
    
    % For next iteration
    t_bar = t_FL_MPC(2:end-1);
    x_bar = x_FL_MPC(:,2:end-1)';
    u_bar = u_lin_k(2:end);
    % TODO Q2) This assumes that the terminal set is the origin.
    x_bar = [x_bar; 0 0];
    u_bar = [u_bar; 0];
    t_bar = [t_bar-dt t_bar(end) + dt];
end

end
