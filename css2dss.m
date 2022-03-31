function [Ad,Bd,Cd] = css2dss(p,Ts,A,B,C)
% Discretize state space model (c2d)
[nx,nu] = size(B);
nd = size(C,2);

switch p.Disc
    case 'Exact'
        dss = expm([A B C; zeros(nu+nd,nx+nu+nd)]*Ts);
        Ad = dss(1:nx,1:nx);
        Bd = dss(1:nx,nx+1:nx+nu);
        Cd = dss(1:nx,nx+nu+1:nx+nu+nd);
    case 'Forward Euler'
        Ad = eye(size(A))+A*Ts;
        Bd = B*Ts;
        Cd = C*Ts;
    case 'Backward Euler'
        Ad = inv(eye(size(A))-A*Ts);
        Bd = Ad*B*Ts;
        if ~isempty(C)
            Cd = Ad*C*Ts;
        else
            Cd = [];
        end
    case 'Tustin'
        Ad = (eye(size(A))+A*Ts/2)*inv(eye(size(A))-A*Ts/2);
        Bd = inv(eye(size(A))-A*Ts/2)*B*(Ts);
        if ~isempty(C)
            Cd = inv(eye(size(A))-A*Ts/2)*C*(Ts);
        else
            Cd = [];
        end
end
end