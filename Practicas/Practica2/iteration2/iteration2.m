%%Segunda Iteracion
%5.1 Aumentamos los rangos
T.Range('theta_h',1) = 0.5;
T.Range('lambda',1) = 1500;
T.Range('mu_m',1) = 0.12;
T.Range('Ms0',2) = 1.5e6;
T.Range('theta_m',2) = 1.2;
save('iteration1/values/T.mat','T')

%%
%5.2 estimation
opt = optimoptions('lsqcurvefit','UseParallel',true,'Display','iter');
[T7_1,res_1] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); 
save('iteration2/values/Results7.mat','T7_1','res_1','xdata','ydata2')

%%
%5.2 Identifiability analysis 
%Como hay solo una familia se cogen todas las estimaciones
th_1 = sum(res_1<res_1(1)*1.01)
y7_1 = gsua_eval(T7_1.Estlsqc(:,1:th_1),T7_1,xdata,ydata2);
savefig('iteration2/figures/curves.fig')

%%
bestest_1 = y7_1(1,:);
trend = [bestest_1(1),bestest_1(2:end)-bestest_1(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration2/figures/EstimatedvsReal.fig')

%%
T7_1.Nominal = T7_1.Estlsqc(:,1);
%res : funciones de costo

T7_1ia = gsua_ia(T7_1,T7_1.Estlsqc(:,1:th_1), false, true); 
%savefig('iteration2/figures/CorrelationsDiag.fig')

