function [P_trace, K_norm, x_est] = KalmanExtended(A,C,odes,t,x,u,y,x0,M,R,P0,Ts,V)
%KALMAN Summary of this function goes here
%   Detailed explanation goes here

%Inicializando variables importantes
n = length(y); x_est = zeros(n,2); P_trace = zeros(n,1); K_norm = zeros(n,1);

if nargin < 13
    V = ones(n,2);
end
%Definiendo A como una funcion para poder evaluar facilmente en ella
syms T
A = eye(2) + A*T; A = matlabFunction(A);
%A(T,x1,x2): función de la matriz del modelo linealizado y discretizado

x_act = x0; P_act = P0;

try
    C = double(C);
    C_i = C;
    functionC = false;
catch
    C = matlabFunction(C);
    functionC = true;
end

for k=1:n
    if functionC
        C_i = C(x_act(1),x_act(2));
    end
    A_i = A(Ts,x_act(1),x_act(2));

    % Predicción:
    try
        x_pred = x_act + Ts*odes(u(k),x_act(1),x_act(2));
    catch
        x_pred = x_act + Ts*odes(x_act(1),x_act(2));
    end

    P_pred = A_i*P_act*A_i' + V(k,:)*M*V(k,:)';
    % Corrección:
    K = (P_pred*C_i')*inv(C_i*P_pred*C_i' + R);
    x_act = x_pred + K*(y(k) - C_i*x_pred);
    P_act = (eye(2) - K*C_i)*P_pred;
    % Estimación:
    x_est(k,:) = x_act';
    % Medidas:p
    P_trace(k) = trace(P_pred);
    K_norm(k) = norm(K);
end

end

