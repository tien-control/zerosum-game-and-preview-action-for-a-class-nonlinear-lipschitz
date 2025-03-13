%% Nonzero sum games
clc; clear; close all;
%%  The Flexible-link robot model
A = [0 1 0 1
    -48.6 -1.25 48.6 0
    0 0 0 1
    1.95 0 -1.95 0];
B = [0; 21.6; 0; 0];
C = [0 0 1 0];
E = [1; 1; 1; 1];

n = size(A,2);
m = size(B,2);
p = size(C,1);
q = size(E,2);

Q_e = 0*eye(p);
Q_x = 1*eye(n);
R = 1*eye(m);

%% Time step
Step = 0.001;
T_end = 80;
t = 0:Step:T_end;

%% Output tracking
A_n0 = [zeros(q,p) C; zeros(n,q) A];
B_n = [zeros(p,m); B];

E_n = [zeros(p,q) -eye(p); E zeros(n,p)];
M_n = zeros(n+q);
N_n = zeros(n+q);

%% Parameter
gamma = 8;
a1 = 1;
a2 = 1;
a3 = 1;
F1 = 1*eye(10);
F2 = 1*eye(10);
F3 = 1*eye(10);
F4 = 1*eye(10);
Q = blkdiag(Q_e,Q_x);
xx = [];
WW = [];
%% Initial value
x = 0.1*rand(5,1);
W1 = 0.1*rand(10,1);
W2 = 0.1*rand(10,1);
W3 = [0;0;0;0;0;0;0;0;0;0];

%% Simulation
    for i = 1:length(t)
        t(i)
        noise = 0.1*(sin(5*pi*t(i)) - sin(exp(1)*t(i)) + sin(t(i))^5 - cos(20*t(i))^5 + sin(-1.2*t(i))^2*cos(0.5*t(i)));

        u_m = -1/2*(R^-1)*B_n'*d_phi(x)'*W2;
        d_m = 1/(2*gamma^2)*E_n'*d_phi(x)'*W3;

        N_n(p+1:end,p+1:end) = df_fun(x);
        A_n = A_n0 + N_n;

        dx = A_n*x + B_n*(u_m+noise) + E_n*(d_m+noise);

        sigma = d_phi(x)*(A_n*x + B_n*u_m + E_n*d_m);
        dW1 = -a1*sigma/((sigma'*sigma+1)^2)*(sigma'*W1+x'*Q*x-gamma^2*norm(d_m)^2+u_m'*R*u_m);

        Da1 = d_phi(x)*B_n*(R^-1)*B_n'*d_phi(x)';
        Ea1 = d_phi(x)*E_n*E_n'*d_phi(x)';
        m = sigma/((sigma'*sigma+1)^2);        
        sigma_n = sigma/(sigma'*sigma+1);
        
        V_t = W1'*phi(x);
        W1 = W1 + Step*dW1;
        V_s = W1'*phi(x);


        dW2 = -a2*((F2*W2-F1*W1)-(1/4)*Da1*W2*m'*W1);
        dW3 = -a3*((F4*W3-F3*W1)+1/(4*gamma^2)*Ea1*W3*m'*W1);

        W2 = W2 + Step*dW2; 
        W3 = W3 + Step*dW3;
        
        if norm(V_s - V_t) < 1e-12 && t(i) > 10
            "Break"
            t(i)
            break
        end
        
        if i == length(t)
            break
        end
        %% Update state
        x = x + Step*dx;

        xx = [xx,x];
        WW = [WW,W1];
    end

    for j = i:i+50000
        u_m = -1/2*(R^-1)*B_n'*d_phi(x)'*W2;
        d_m = 1/(2*gamma^2)*E_n'*d_phi(x)'*W3;
    
        N_n(p+1:end,p+1:end) = df_fun(x);
        A_n = A_n0 + N_n;
        dx = A_n*x + B_n*u_m + E_n*d_m;
    
        %W1 = W1 + Step*dW1;
        %W2 = W2 + Step*dW2;
        %W3 = W3 + Step*dW3;
            %% Update state
        x = x + Step*dx;
        xx = [xx,x];
        WW = [WW,W1];
    end

figure(1);
plot(1:size(xx,2),xx);
legend('$de$','$d^2x_1$','$d^2x_2$','$d^2x_3$','$d^2x_4$', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 10);

figure(2);
plot(1:size(WW,2),WW);

W1
W2
W3

%% Cac ham khac

function a = phi(x)
a = [x(1)^2; x(2)^2; x(3)^2; x(4)^2; x(5)^2; x(1)^4; x(2)^4; x(3)^4; x(4)^4; x(5)^2];
end

function a = d_phi(x)
a = [2*x(1)  0   0   0   0
     0   2*x(2)  0   0   0
     0   0   2*x(3)  0   0
     0   0   0   2*x(4)  0
     0   0   0   0   2*x(5)
     4*x(1)^3    0   0   0   0
     0   4*x(2)^3    0   0   0
     0   0   4*x(3)^3    0   0
     0   0   0   4*x(4)^3    0
     0   0   0   0   4*x(5)^3 ];
end         

%% Ham bo sung
function a = f_fun(x)
a = [0
     0
     0
    -0.33*sin(x(3))];
end

function a = df_fun(x)
a = [0 0 0 0
     0 0 0 0
     0 0 0 0
     0 0 -0.33*cos(x(3)) 0];
end