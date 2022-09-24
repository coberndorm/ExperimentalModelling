clc
clear

out = sim('Dos_tanques_acoplados.slx');
yexp = out.h2;
uexp = out.u;
texp = out.tout;

%% Parámetros
e = 1.1; % error de modelado
a = 0.0005*e; K1 = 0.047*e; K2 = 0.065*e; A1 = 0.45*e; %A2 = 0.45;

T = 5;

%% Condiciones iniciales
I = eye(3);
x_pre = [1.2; 0.8; 1];
P_pre = 0.3*I;
Q = 0.2*I;
R = 1;
G = 1;

%% 
x_est = []; P = []; K0 = [];
C = [0 1 0];
N = 1;
n = length(yexp);

% Monitorear número condición y rango de la matriz de observabilidad
condi = [];
ranks = [];

% Intervalos
x_est_sup = [];
x_est_inf = [];

for k = 0:n-1
    
    % Corrección cada N 
    if k/N == round(k/N)
        K = P_pre*C'*inv(C*P_pre*C' + R);
        x_act = x_pre + K*(yexp(k+1) - C*x_pre);
        P_act = (I - K*C)*P_pre;
    else
        P_act = P_pre;
        x_act = x_pre;
    end
    % Predicción
    x1 = x_act(1);
    x2 = x_act(2);
    x3 = x_act(3);

    x_pre = [x1 - (K1*T/A1)*sqrt(x1-x2) + a*T*uexp(k+1)/A1;
                x2 + K1*T*sqrt(x1-x2)/x3 - K2*T*sqrt(x2)/x3;
                x3];

    Fi = [1 - (K1*T)/(2*A1*sqrt(x1-x2)), K1*T/(2*A1*sqrt(x1-x2)), 0; 
        -K1*T/(2*x3*sqrt(x1-x2)), 1 + K1*T/(2*x3*sqrt(x1-x2)) - K2*T/(2*x3*sqrt(x2)), ...
        -K1*T*sqrt(x1-x2)/(x3^2) + K2*T*sqrt(x2)/(x3^2); 0,0,1];

    P_pre = Fi*P_act*Fi' + G*Q*G';
    x_est = [x_est x_pre];
    P = [P trace(P_pre)];
    K0 = [K0 norm(K)];

    M0 = [C; C*Fi; C*(Fi^2)];
    condi = [condi cond(M0)];
    ranks = [ranks rank(M0)];

    x_est_sup = [x_est_sup x_pre + 1.96*sqrt(diag(P_pre))];
    x_est_inf = [x_est_inf x_pre - 1.96*sqrt(diag(P_pre))];
end

%% Plot
t = T*[0:(n-1)];

figure(1), plot(t,yexp,t,x_est(2,:)),legend({'exp','estimado'}), title('x2')
hold on, plot(t, x_est_sup(2,:), t, x_est_inf(2,:)), 
legend({'exp','estimado','sup','inf'}), hold off

figure(2), plot(t,x_est(1,:)), title('x1 estimado')
hold on, plot(t, x_est_sup(1,:), t, x_est_inf(1,:)), 
legend({'estimado','sup','inf'}), hold off

figure(3), plot(t,P), title('Traza de P')
figure(4), plot(t,K0), title('Norma de K')
figure(5), plot(t,condi), title('Número condición')
figure(6), plot(t,ranks), title('Rango Matriz de Observabilidad')

figure(7), plot(t,x_est(3,:), t, x_est_sup(3,:),t,x_est_inf(3,:)), 
%figure(7), plot(t,x_est(3,:)),
title('Estimación parámetro x3=A2')

%% Error cuadrático medio del estado estimado versus experimental
m = sum((x_est(2,:)-yexp').^2)/n
%% Residuos
figure(8), plot(t, yexp' - C*x_est), title('Residuos x2')