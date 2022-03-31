%% Linearize about the origin

dt = 1;

p = params();
[dyn, o] = dynamics_and_outputs(p);

A_origin = dyn.Df_func_model(0,0)+dyn.Dg_func_model(0,0)*0;
B_origin = dyn.g_func_model(0,0);
C_origin = dyn.f_func_model(0,0) + dyn.g_func_model(0,0)*0 - (A_origin*[0;0] + B_origin*0);
[Ad_,Bd_,Cd_] = css2dss(p,dt,A_origin,B_origin,C_origin);
Ad_at_origin = Ad_;
Bd_at_origin = Bd_;
Cd_at_origin = Cd_;

%% Ugo's code
% A = Ad_at_origin;
% B = Bd_at_origin;
% Q = eye(2);
% R = 1;
% 
% X = [100 100 -100 -100;100 -100 100 -100];
% U = [-100 100];
% 
% [P,L,G] = dare(A, B, Q,R); 
% Finf=-(B'*P*B+R)^(-1)*B'*P*A;
% Acl = A+B*Finf;
% Ptr = dlyap(Acl',Finf'*R*Finf+Q);
% system = LTISystem('A',Acl);
% Xtilde = Polyhedron('A',[eye(2);-eye(2);Finf;-Finf],'b',[X.V(1,1);X.V(1,1);X.V(1,1);X.V(1,1);U.V(1,1);U.V(1,1)]);
% Oinf = system.invariantSet('X',Xtilde);
% Af = Oinf.A;
% bf = Oinf.b;


%% Solve LQR to find a feedback gain

K = lqrd(Ad_at_origin, Bd_at_origin, eye(2), 1, [0; 0], dt);

A_cl = Ad_at_origin - Bd_at_origin*K;


%% Calculate Backward Reachable Set

F = [1 0; 0 1; -1 0; 0 -1];
b = [0.1 0.1 0.1 0.1]';

F = F*A_cl;

[X,Y] = meshgrid(linspace(-10,10));
V = zeros(size(X));
for i=1:numel(X)
    V(i) = all(F*[X(i); Y(i)]-b <= 0);
end

surf(V)


%% Verify feasibility of vertices