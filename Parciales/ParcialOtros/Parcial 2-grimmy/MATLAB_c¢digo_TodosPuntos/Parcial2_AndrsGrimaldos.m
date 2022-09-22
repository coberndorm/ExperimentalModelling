% Andrés Grimaldos Echavarría
%% Punto 1
Datos = [0 1 0.0052;
         1 1 2.1190;
         2 1 1.4082;
         3 1 1.1868;
         4 1 1.2544;
         5 1 1.2402;
         6 1 1.2463;
         7 1 1.2448];
[phi, Y, param, cova] = MCO(Datos, 1, 4, 1);
% Debido a la forma de phi'*phi, nos quedan columnas iguales, luego no es
% invertible, por tanto, debemos modificar las columnas que usan elementos
% de las entradas, usando instantes pasados de tal manera que se conviertan
% en 0 y para un arx121, de la siguiente manera:
[phi, Y, param, cova] = MCO(Datos, 1, 2, 1);
phi(1,3) = 0;
param = inv(phi'*phi)*(phi')*Y;
data = iddata(Y, phi(:,end),1);
orden_excitacion = pexcit(data);
V = 0.5*(Y'*Y - Y'*phi*inv(phi'*phi)*phi'*Y);
lambda = 2*V/(8-1-3-1);
cova = lambda*inv(phi'*phi);
desv = sqrt(diag(cova));
G = tf([param(2) param(3)], [1 param(1)], 'Ts', 1, 'InputDelay', 1);
step(G, 0:1:7);
hold on
stairs(0:1:7, Datos(:,3),'r');
legend({'Estimación', 'Experimentales dados'}, 'FontSize', 15, 'Location', 'best');
xlabel('$t$', 'Interpreter', 'Latex', 'FontSize', 15); 
ylabel('$y(t)$', 'Interpreter', 'Latex', 'FontSize', 15);
title('Comparación FDT estimada con datos experimentales dados', 'FontSize', 10);
grid on 
%% Punto 2
syms x1 x2 u a b c d 
vars = [x1, x2];
odes = [-a*x1 + b*x1*x2; c*x1^2 - d*x2 + u];
% Revisamos que el Jacobiano esté bien:
J = jacobian(odes, vars);
syms T
A = eye(2) + J*T;
%A(T,x1,x2): función de la matriz del modelo linealizado y discretizado
A = matlabFunction(A); 
Ts = 0.1; % tiempo de muestreo T.
tiempo = 0:Ts:50; % periodo de tiempo a usar.
x0 = [1, 1]; % Condición inicial
% Probamos dos entradas:
u = sin(5*tiempo) + randn(1,length(tiempo));%ones(1, length(tiempo));%sin(5*tiempo) + randn(1,length(tiempo));
plot(tiempo, u)
title('Entrada');
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$u(k)$', 'Interpreter', 'Latex', 'FontSize', 15);
grid on
a = 0.01; b = 0.02; , c = 0.1; d = 0.2;
[t, x] = ode45(@(t,x) ode_parcial(t, x, u, tiempo, a, b, c, d), tiempo, x0);
y = x(:,1);
% Condiciones iniciales dados por el usuario:
M = 0.01*eye(2); N = 0.5; P0 = M; x0 = [1 1]'; Ts = 0.1;
% Parámetros conocidos:
C = [1 0];
n = length(y); x_est = zeros(n,2); P_trace = zeros(n,1); K_norm = zeros(n,1);
x_act = x0; P_act = P0;
for k = 1:n
    A_ = A(Ts,a,b,c,d,x_act(1), x_act(2));
    % Predicción:
    x_pred = [x_act(1) + Ts*(b*x_act(1)*x_act(2)-a*x_act(1));
              x_act(2) + Ts*(c*x_act(1)^2 - d*x_act(2) + u(k))];
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
plot(t, x(:,2), t, x_est(:,2))
legend({'Modelo original', 'Estimación Kalman Extendido'}, 'Location', 'best');
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$x_{2}(k)$', 'Interpreter', 'Latex', 'FontSize', 15);
title('Comparación modelo original con estimación kalman extendido', 'FontSize', 12);
grid on