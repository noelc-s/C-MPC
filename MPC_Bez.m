function [X_Lin_MPC, T_Lin_MPC, U_Lin_MPC, X_K_MPC_CLF, XD_Lin_MPC, u_Lin_MPC, U_FF_MPC_CLF] = MPC_Bez(p, dyn, o, T_FL_MPC, X_FL_MPC, U_FL_MPC)
disp('Setting up proposed approach');

% TODO: Q1) need to use a feasible (x,u) to start, which the following is not.
t_bar = T_FL_MPC;
x_bar = X_FL_MPC;
u_bar = U_FL_MPC;

x = p.ODE.X0;
X_Lin_MPC = [];
XD_Lin_MPC = [];
T_Lin_MPC = [];
U_Lin_MPC = [];
X_K_MPC_CLF = [];
U_FF_MPC_CLF = [];

N = p.MPC.N;
Xf = p.ODE.Xf;
dt = p.MPC.dt;
low_level_dt = p.ll_dt;
alpha_MPCFL = p.MPC_CLF.alpha_MPCFL; % Lf
beta_MPCFL = p.MPC_CLF.beta_MPCFL;  % Lg
gamma_MPCFL = p.MPC_CLF.gamma_MPCFL; % Robust_level_Set_value
Gamma_MPCFL = p.MPC_CLF.Gamma_MPCFL; % E_bar (i.e. sqrt(gamma/lambda_min(P))
delta_MPCFL = p.MPC_CLF.delta_MPCFL;
state_stage_cost = p.MPC_CLF.state_stage_cost;
input_stage_cost = p.MPC_CLF.input_stage_cost;
aux_stage_cost = p.MPC_CLF.aux_stage_cost;
A_in = p.Const.A_in;
b_in = p.Const.b_in;
P_lyap = p.CLF.P_lyap;

order = p.MPC_CLF.order;
HP_lyap = [diag(-ones(order,1)) zeros(order,1)] + [zeros(order,1) diag(ones(order,1))];
HP_lyap = HP_lyap*order;
HP_lyap2 = HP_lyap(1,:);
for i = 1:order-1
    HP_lyap2 = [HP_lyap2;
        i/order*HP_lyap(i,:) + (order-i)/order*HP_lyap(i+1,:)];
end
HP_lyap2 = [HP_lyap2; HP_lyap(order,:)];
H = HP_lyap2/dt;
H2 = H^2;

M = [2*alpha_MPCFL*beta_MPCFL beta_MPCFL;beta_MPCFL 0];
[evec, eval] = eig(M);
% normalized_evec = evec(:,1)/norm(evec(:,1));
M_k = eval(2,2)*evec(:,2)*evec(:,2)';
if(any(eig(M_k) < -eps))
   error('Construction of M_k failed'); 
end

% M_k = beta_MPCFL*[delta_MPCFL^3 delta_MPCFL^2;
%     delta_MPCFL^2 delta_MPCFL];

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
N_k = zeros(2,N-1);
Gamma_k = zeros(N-1,1);
x_lin_k = zeros(N-1,2);
u_lin_k = zeros(N-1,1);
Ad_km1 = zeros(state_dim, state_size_without_terminal);
Bd_km1 = zeros(state_dim, N-1);
Cd_km1 = zeros(state_dim, N-1);
N_km1 = zeros(2,N-1);
Gamma_km1 = zeros(N-1,1);
x_lin_km1 = zeros(N-1,2);
u_lin_km1 = zeros(N-1,1);

% Optimizer variables
states = sdpvar(state_size_with_terminal,1);
inputs = sdpvar(input_size,1);
auxiliaries = sdpvar(auxiliary_size,1);
slack = sdpvar(1);
Ad_ = sdpvar(state_dim, state_size_without_terminal);
Bd_ = sdpvar(state_dim, N-1);
Cd_ = sdpvar(state_dim, N-1);
N_k_ = sdpvar(2,N-1);
Gamma_k_ = sdpvar(N-1,1);
x_lin_ = sdpvar(N-1,2);
u_lin_ = sdpvar(N-1,1);
x0_ = sdpvar(state_dim,1);

% adding a slack variable for state constraints
Constraints = [];

% clear pos1 pos2 pos3 pos4 vel1 vel2 vel3 vel4 u_k s
pos1 = sdpvar(N-1,1);
pos2 = sdpvar(N-1,1);
pos3 = sdpvar(N-1,1);
pos4 = sdpvar(N-1,1);

vel1 = sdpvar(N-1,1);
vel2 = sdpvar(N-1,1);
vel3 = sdpvar(N-1,1);
vel4 = sdpvar(N-1,1);

u_k = sdpvar(N-1,1);
s = sdpvar(2,N-1);
t = sdpvar(1,N-1);
w = sdpvar(3,N-1);

[Low, Diag] = ldl(M_k);
P_ldl = real(sqrtm(Diag))*Low'; % it's real, I swear

for i = 1:N-1
    % Decision variables
    Constraints = [Constraints pos1(i,1) == states((i-1)*state_dim+1)];
    Constraints = [Constraints vel1(i,1) == states((i-1)*state_dim+2)];
    Constraints = [Constraints pos4(i,1) == states((i)*state_dim+1)];
    Constraints = [Constraints vel4(i,1) == states((i)*state_dim+2)];
    % Defined:
    Constraints = [Constraints pos2(i,1) == dt/3*vel1(i,1)+pos1(i,1)];
    Constraints = [Constraints pos3(i,1) == -dt/3*vel4(i,1)+pos4(i,1)];
    Constraints = [Constraints vel2(i,1) == H(2,:)*[pos1(i,1) pos2(i,1) pos3(i,1) pos4(i,1)]'];
    Constraints = [Constraints vel3(i,1) == H(3,:)*[pos1(i,1) pos2(i,1) pos3(i,1) pos4(i,1)]'];
    
    Constraints = [Constraints u_k(i,1) == inputs((i-1)*input_dim+1:i*input_dim)];
    Constraints = [Constraints s(:,i) == auxiliaries((i-1)*auxiliary_dim+1:(i)*auxiliary_dim)];
    
    % Dynamics
    Constraints = [Constraints [pos4(i,1); vel4(i,1)] == Ad_(1:state_dim,(i-1)*state_dim+1:(i)*state_dim)*[pos1(i,1); vel1(i,1)] ...
        + Bd_(1:state_dim, i)*u_k(i,1) + Cd_(1:state_dim, i)];
    %     Constraints = [Constraints norm(sqrtm(M_k)*s(:,i) + 1/2*inv(sqrtm(M_k))*N_k_(:,i),2) <= (1/4*N_k_(:,i)'*inv(M_k)*N_k_(:,i)-(u_max- Gamma_k_(i)))^(1/2)];
    %     Constraints = [Constraints s(:,i)'*M_k*s(:,i) + N_k_(:,i)'*s(:,i) + Gamma_k_(i) <= u_max];
    Constraints = [Constraints w(:,i) == [s(:,i); t(:,i)]];
    Constraints = [Constraints t(:,i) + 1/4 + N_k_(:,i)'*s(:,i) + Gamma_k_(i) <= p.Const.u_max];
    Constraints = [Constraints norm([P_ldl zeros(2,1); zeros(1,2) 1]*w(:,i),2) <= t(:,i) + 1/2];
    Constraints = [Constraints [norm([pos1(i,1); vel1(i,1)]-x_lin_(i,:)',2); ...
        norm([vel1(i,1) vel2(i,1) vel3(i,1) vel4(i,1)]*H(1,:)'-o.Lf2y_func(x_lin_(i,1),x_lin_(i,2)),2)] <= s(:,i)];
    Constraints = [Constraints [norm([pos2(i,1); vel2(i,1)]-x_lin_(i,:)',2); ...
        norm([vel1(i,1) vel2(i,1) vel3(i,1) vel4(i,1)]*H(2,:)'-o.Lf2y_func(x_lin_(i,1),x_lin_(i,2)),2)] <= s(:,i)];
    Constraints = [Constraints [norm([pos3(i,1); vel3(i,1)]-x_lin_(i,:)',2); ...
        norm([vel1(i,1) vel2(i,1) vel3(i,1) vel4(i,1)]*H(3,:)'-o.Lf2y_func(x_lin_(i,1),x_lin_(i,2)),2)] <= s(:,i)];
%     Constraints = [Constraints inputs((i-1)*input_dim+1:i*input_dim) <= p.Const.u_max-Gamma_k_(i)];
%     Constraints = [Constraints inputs((i-1)*input_dim+1:i*input_dim) >= p.Const.u_min+Gamma_k_(i)];
    % b_k in X \ominus E
    Constraints = [Constraints A_in*[pos1(i,1); vel1(i,1)] <= b_in-sqrt(2*gamma_MPCFL*p.noise_mag^2)*diag(A_in*P_lyap*A_in')];%+slack];
    Constraints = [Constraints A_in*[pos2(i,1); vel2(i,1)] <= b_in-sqrt(2*gamma_MPCFL*p.noise_mag^2)*diag(A_in*P_lyap*A_in')];%+slack];
    Constraints = [Constraints A_in*[pos3(i,1); vel3(i,1)] <= b_in-sqrt(2*gamma_MPCFL*p.noise_mag^2)*diag(A_in*P_lyap*A_in')];%+slack];
    Constraints = [Constraints A_in*[pos4(i,1); vel4(i,1)] <= b_in-sqrt(2*gamma_MPCFL*p.noise_mag^2)*diag(A_in*P_lyap*A_in')];%+slack];
    Constraints = [Constraints s(:,i) >= 0];
    %     Constraints = [Constraints norm(s(:,i),2) <= 1];
end

% Terminal Set
Constraints = [Constraints states(state_dim*(N-1)+1:state_dim*(N-1)+2) == Xf];

% Add this for ISS tube, which is not considered yet.
Constraints = [Constraints V_point(states(1:state_dim)-x0_) <= gamma_MPCFL*p.noise_mag^2];
% Constraints = [Constraints states(1:2) == x0_]; % placeholder

Q = state_stage_cost*eye(state_dim*N);
R = input_stage_cost*eye(input_dim*(N-1));
Aux = aux_stage_cost*eye(auxiliary_dim*(N-1));
W = 1e9;
% Q = blkdiag(Q,R,Aux,1e9);
% f = zeros(total_dim+1,1);
Objective = 1/2*(states'*Q*states + inputs'*R*inputs + auxiliaries'*Aux*auxiliaries + slack'*W*slack);

% SP_lyapQ_epsilon = 0.1;
% SQP_lyapconverged = false;

P = optimizer(Constraints,Objective,sdpsettings('solver','mosek','verbose',0),...
    {Ad_,Bd_,Cd_,N_k_,Gamma_k_,x_lin_,u_lin_,x0_},...
    {states,inputs,auxiliaries,slack,pos1,pos2,pos3,pos4,vel1,vel2,vel3,vel4,u_k,s,t,w});

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
norm_G_pinv = norm(pinv(o.LgLfy_func(x_F(i,1),x_F(i,2))),2);
norm_K = norm(p.CLF.alpha_CLF,2);
N_k_at_origin = [2*alpha_MPCFL*beta_MPCFL*Gamma_MPCFL+alpha_MPCFL*norm_G_pinv+beta_MPCFL*norm_K*Gamma_MPCFL, ...
            norm_G_pinv + beta_MPCFL*Gamma_MPCFL];
Gamma_k_at_origin = (beta_MPCFL*Gamma_MPCFL^2+norm_G_pinv*Gamma_MPCFL)*(alpha_MPCFL + norm_K);

delay_buffer = repmat(x(end,:)', 1, p.ll_delay);

for iter = 1:p.ODE.tspan(end)/dt
    disp(iter)
    sampled_x0 = delay_buffer(:,1);
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
        norm_G_pinv = norm(pinv(o.LgLfy_func(x_lin_k(i,1),x_lin_k(i,2))),2);
        norm_K = norm(p.CLF.alpha_CLF,2);
        N_k(:,i) = [2*alpha_MPCFL*beta_MPCFL*Gamma_MPCFL+alpha_MPCFL*norm_G_pinv+beta_MPCFL*norm_K*Gamma_MPCFL, ...
            norm_G_pinv + beta_MPCFL*Gamma_MPCFL];
        Gamma_k(i) = (beta_MPCFL*Gamma_MPCFL^2+norm_G_pinv*Gamma_MPCFL)*(alpha_MPCFL + norm_K);
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
    [sol, diagnostics,d1,d2,d3,d4] = P({Ad_k,Bd_k,Cd_k,N_k,Gamma_k,x_lin_k,u_lin_k,sampled_x0});
    if iter == 1 & diagnostics ~= 0
        error('Issue with Mosek in proposed');
    elseif diagnostics ~= 0
        warning('Falling back on previous linearization')
        Ad_km1 = [Ad_km1(:,3:end) Ad_at_origin];
        Bd_km1 = [Bd_km1(:,2:end) Bd_at_origin];
        Cd_km1 = [Cd_km1(:,2:end) Cd_at_origin];
        N_km1 = [N_km1(:,2:end) N_k_at_origin'];
        Gamma_km1 = [Gamma_km1(2:end); Gamma_k_at_origin];
        % TODO: assumes terminal point is unforced equilibrium
        x_lin_km1 = [x_lin_km1(2:end,:); 0 0];
        u_lin_km1 = [u_lin_km1(2:end,:); 0];
        [sol, diagnostics,d1,d2,d3,d4] = P({Ad_km1,Bd_km1,Cd_km1,N_km1,Gamma_km1,x_lin_km1,u_lin_km1,sampled_x0});
    else
        Ad_km1 = Ad_k;
        Bd_km1 = Bd_k;
        Cd_km1 = Cd_k;
        N_km1 = N_k;
        Gamma_km1 = Gamma_k;
        x_lin_km1 = x_lin_k;
        u_lin_km1 = u_lin_k;
    end
    
    %     MPC_sol = value(vars);
    
    t_FL_MPC = 0:dt:dt*(N-1);
    %     x_FL_MPC = [MPC_sol(1:2:state_end_index) MPC_sol(2:2:state_end_index)];
    %     u_lin = MPC_sol(state_end_index+1:input_end_index);
    %     aux_vars = MPC_sol(input_end_index+1:auxiliary_end_index);
    %     slack = MPC_sol(end);
    x_FL_MPC = [sol{1}(1:2:state_end_index) sol{1}(2:2:state_end_index)];
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
    a = x_FL_MPC(1,1);
    b = 1/3*(dt*x_FL_MPC(1,2)+3*x_FL_MPC(1,1));
    c = 1/3*(-dt*x_FL_MPC(2,2)+3*x_FL_MPC(2,1));
    d = x_FL_MPC(2,1);
    d_x = @(t) a*(t/dt).^0.*(1-t/dt).^3+3*b*(t/dt).^1.*(1-t/dt).^2+3*c*(t/dt).^2.*(1-t/dt).^1+d*(t/dt).^3.*(1-t/dt).^0;
    d_x_d = @(t) (3*d*t.^2)/dt.^3 - (3*c*t.^2)/dt.^3 - (3*a*(t/dt - 1).^2)/dt + (3*b*(t/dt - 1).^2)/dt + (6*b*t.*(t/dt - 1))/dt^2 - (6*c*t.*(t/dt - 1))/dt^2;
    d_x_dd = @(t) (12*b*(t/dt - 1))/dt^2 - (6*a*(t/dt - 1))/dt^2 - (6*c*(t/dt - 1))/dt^2 + (6*b*t)/dt^3 - (12*c*t)/dt^3 + (6*d*t)/dt^3;
    
    %%%%%%%%%%%%%%%% Which low level controller to use %%%%%%%%%%%%%%%%%%%
    F = [0 1; 0 0];
    G = [0; 1];
    A_cl = F - G*p.CLF.alpha_CLF;
    Q_lyap = eye(2);
    P_lyap_lyap = lyap(A_cl', Q_lyap);
    gamma = 1/max(eig(P_lyap_lyap));
    eta = @(x,t) [x(1) - d_x(t); x(2) - d_x_d(t)];
    V = @(x,t) eta(x,t)'*P_lyap_lyap*eta(x,t);
    LFV = @(x,t) eta(x,t)'*(F'*P_lyap_lyap + P_lyap_lyap*F)*eta(x,t);
    LGV = @(x,t) 2*eta(x,t)'*P_lyap_lyap*G;
    options = optimset('Display', 'off');
    sigma = 1;
    
    %%% FL CLF1
    %     v = @(x,t) quadprog(1,0, LGV(x,t), -gamma*V(x,t)-LFV(x,t),[],[],[],[],[],options);
    v = @(x,t) -max(0,(LFV(x,t)+gamma*V(x,t))/(LGV(x,t)'*LGV(x,t)))*LGV(x,t)';
    FL_CLF1 = @(t,x) o.LgLfy_func(x(1),x(2))\(-o.Lf2y_func(x(1),x(2)) + v(x,t) + d_x_dd(t));
    
    %%% FL CLF2
    calF = @(x,t) Lf2y_func(x(1),x(2))-d_x_dd(t);
    calG = @(x,t) LgLfy_func(x(1),x(2));
    v = @(x,t) quadprog(calG(x,t)'*calG(x,t),calF(x,t)'*calG(x,t), LGV(x,t)*calG(x,t), -gamma*V(x,t)-LFV(x,t)-LGV(x,t)*calF(x,t),[],[],[],[],[],options);
    FL_CLF2 = @(t,x) v(x,t);
    
    T = 0;
    X = [];
    for j = 1: dt / low_level_dt
        x0 = x(end,:)';
        switch p.low_level
            case 'None'
                [t,x] = ode45(@(t,x) dyn.f_func_w(x(1),x(2),t) + dyn.g_func_w(x(1),x(2),t)*(u_lin_k(1)),[0 low_level_dt],x0); % no low level
            case 'CLF'
                [t,x] = ode45(@(t,x) dyn.f_func_w(x(1),x(2),t) + dyn.g_func_w(x(1),x(2),t)*(FL_CLF1((j-1)*low_level_dt,x0)),[0 low_level_dt],x0); % CLF1
                %         [t,x] = ode45(@(t,x) f_func(x(1),x(2),t) + g_func(x(1),x(2),t)*(FL_CLF2(t,x')),[0 dt],x0); % CLF2 -- SLOW
            otherwise
                error('That low level controller not implemented yet!');
        end
        T = [T; T(end) + t(2:end)];
        X = [X; x(2:end,:)];
    end
    t = T(2:end);
    x = X;
    delay_buffer(:,1:end-1) = delay_buffer(:,2:end);
    delay_buffer(:,end)= x(end,:)';
    
    
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
        XD_Lin_MPC = [XD_Lin_MPC; [d_x(t) d_x_d(t)]];
        U_Lin_MPC = [U_Lin_MPC; u_Lin_MPC'];
        U_FF_MPC_CLF = [U_FF_MPC_CLF; u_lin_k(1)*ones(size(t))];
    else
        T_Lin_MPC = [T_Lin_MPC; t(2:end)+T_Lin_MPC(end)];
        X_Lin_MPC = [X_Lin_MPC; x(2:end,:)];
        XD_Lin_MPC = [XD_Lin_MPC; [d_x(t(2:end)) d_x_d(t(2:end))]];
        U_Lin_MPC = [U_Lin_MPC; u_Lin_MPC(2:end)'];
        U_FF_MPC_CLF = [U_FF_MPC_CLF; u_lin_k(1)*ones(size(t(2:end)))];
    end
    XD_Lin_MPC(end,:) = NaN; % to separate trajectories
    
    X_K_MPC_CLF = [X_K_MPC_CLF; x_FL_MPC(1,:)];
    
    % For next iteration
    t_bar = t_FL_MPC(2:end-1);
    x_bar = x_FL_MPC(2:end-1,:);
    u_bar = u_lin_k(2:end);
    % TODO Q2) This assumes that the terminal set is the origin.
    x_bar = [x_bar; 0 0];
    u_bar = [u_bar; 0];
    t_bar = [t_bar-dt t_bar(end) + dt];
end

end
