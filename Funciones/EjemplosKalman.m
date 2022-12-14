
%% Ejemplo 2, parcial:
syms x1 x2 u
vars = [x1, x2];
odes = [-0.2*sqrt(x1-x2) + 0.001*u;
        -0.1*sqrt(x2) + 0.2*sqrt(x1-x2)];
%Revisamos el Jacobiano:
J = jacobian(odes, vars);
T = 0.1;
tiempo = 0:T:50;
u = 10*sin(tiempo);%ones(1,length(tiempo)); 
x0 = [3, 2];
[t, x] = ode45(@(t,x) ode2(t, x, u, tiempo), tiempo, x0); 
y = x(:,2);
% Condiciones iniciales dados por el usuario:
M = 0.001*eye(2); N = 0.5; P0 = M; x0 = [3 2]'; Ts = 0.1;
% Parámetros conocidos:
C = [0 1];
n = length(y); x_est = zeros(n,2); P_trace = zeros(n,1); K_norm = zeros(n,1);
x_act = x0; P_act = P0; u = ones(n,1);
for k = 1:n
    A = [1 - 0.1*Ts/sqrt(x_act(1) - x_act(2)) , 0.1/sqrt(x_act(1)-x_act(2));
         0.1/sqrt(x_act(1)-x_act(2)), 1-0.05/sqrt(x_act(2))+T*(0.1/sqrt(x_act(1)-x_act(2)))];
    % Predicción:
    x_pred = [x_act(1) + Ts*(0.001*u(k) - 0.2*sqrt(x_act(1) - x_act(2)));
              x_act(2) + Ts*(-0.1*sqrt(x_act(2))+0.2*sqrt(x_act(1)-x_act(2)))];
    P_pred = A*P_act*A' + M;
    % Corrección:
    K = P_pred*C'*inv(C*P_pred*C' + N);
    x_act = x_pred + K*(y(k) - C*x_pred);
    P_act = (eye(2) - K*C)*P_pred;
    % Estimación:
    x_est(k,:) = x_act';
    % Medidas:
    P_trace(k) = trace(P_act);
    K_norm(k) = norm(K);
end
plot(t, x(:,1), t, x_est(:,1))
legend({'Modelo original', 'Estimación Kalman Extendido'}, 'Location', 'best');
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$x_{1}(k)$', 'Interpreter', 'Latex', 'FontSize', 15);
title('Comparación modelo original con estimación kalman extendido', 'FontSize', 12);
grid on
%% Ejemplo 3, parcial
syms x1 x2 u 
vars = [x1, x2];
odes = [-0.1*x1*x2 + u; 0.1*x1*x2 - 0.9*x2];
% Revisamos que el Jacobiano esté bien:
J = jacobian(odes, vars);
syms T
A = eye(2) + J*T;
A = matlabFunction(A); %A(T,x1,x2)
Ts = 0.1;
tiempo = 0:Ts:50;
x0 = [10, 5];
u = sin(5*tiempo) + randn(1,length(tiempo));%ones(1, length(tiempo));
plot(tiempo, u)
title('Entrada');
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$u(k)$', 'Interpreter', 'Latex', 'FontSize', 15);
grid on
[t, x] = ode45(@(t,x) ode3(t, x, u, tiempo), tiempo, x0);
y = x(:,2);
% Condiciones iniciales dados por el usuario:
M = 0.01*eye(2); N = 0.5; P0 = M; x0 = [10 5]'; Ts = 0.1;
% Parámetros conocidos:
C = [0 1];
n = length(y); x_est = zeros(n,2); P_trace = zeros(n,1); K_norm = zeros(n,1);
x_act = x0; P_act = P0;
for k = 1:n
    A_ = A(Ts, x_act(1), x_act(2));
    % Predicción:
    x_pred = [x_act(1) + Ts*(u(k) - 0.1*x_act(1)*x_act(2));
              x_act(2) + Ts*(0.1*x_act(1)*x_act(2) - 0.9*x_act(2))];
    P_pred = A_*P_act*A_' + M;
    % Corrección:
    K = P_pred*C'*inv(C*P_pred*C' + N);
    x_act = x_pred + K*(y(k) - C*x_pred);
    P_act = (eye(2) - K*C)*P_pred;
    % Estimación:
    x_est(k,:) = x_act';
    % Medidas:
    P_trace(k) = trace(P_act);
    K_norm(k) = norm(K);
end
plot(t, x(:,1), t, x_est(:,1))
legend({'Modelo original', 'Estimación Kalman Extendido'}, 'Location', 'best');
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$x_{1}(k)$', 'Interpreter', 'Latex', 'FontSize', 15);
title('Comparación modelo original con estimación kalman extendido', 'FontSize', 12);
grid on