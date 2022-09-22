%% Corrección punto 1:
Datos = [0 1 0.0052;
         1 1 2.1190;
         2 1 1.4082;
         3 1 1.1868;
         4 1 1.2544;
         5 1 1.2402;
         6 1 1.2463;
         7 1 1.2448];
[phi, Y, param, cova] = MCO(Datos, 1, 4, 1);
R = phi'*phi/(9);
determ = det(R);
phi = zeros(7,4);
phi(:,1) = Datos(1:end-1,3);
phi(:,2) = ones(7,1);
phi(2:end,3) = ones(6,1);
phi(3:end, 4) = ones(5,1);
phi(4:end, 5) = ones(4,1);
Y = Datos(2:end,3);
data = iddata(Y,phi(:,end),1);
orden_excit = pexcit(data);
param = inv(phi'*phi)*(phi')*Y;
V = 0.5*(Y'*Y - Y'*phi*inv(phi'*phi)*phi'*Y);
lambda = 2*V/(7-4-1);
cova = lambda*inv(phi'*phi);
desv = sqrt(diag(cova));
G = tf([param(2) param(3) param(4) param(5)], [1 param(1) 0 0 0], 'Ts', 1);
step(G, 0:1:7);
hold on
stairs(0:1:7, Datos(:,3),'r');
legend({'Estimación', 'Experimentales dados'}, 'FontSize', 15, 'Location', 'best');
xlabel('$t$', 'Interpreter', 'Latex', 'FontSize', 15); 
ylabel('$y(t)$', 'Interpreter', 'Latex', 'FontSize', 15);
title('Comparación FDT estimada con datos experimentales dados', 'FontSize', 10);
grid on 