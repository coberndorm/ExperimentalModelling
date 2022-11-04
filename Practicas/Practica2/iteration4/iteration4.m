%%Cuarta Iteracion
%Fijamos Beta_h


Range(9,:) = T7_3.Nominal('beta_h');
[T,~]=gsua_dpmat(odes,vars,[0 150],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output = 1;
M = gsua_dmatrix(T,100);

%Realizamos analisis otra vez
Tsa = gsua_sa(M,T,'parallel', false, 'SensMethod', 'Xiao', 'ynom', ynom);

gsua_plot('Bar',Tsa,Tsa.STi)
savefig('iteration4/figures/SensibilityAnalisis.fig')


% Confiabilidad

c = sum(Tsa.Si)/sum(abs(Tsa.Si))

%5.2 estimation
opt = optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7_3,res_3] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); 
save('iteration4/values/Results7.mat','T7_3','res_3','xdata','ydata2')

%5.2 Identifiability analysis 
th_3 = sum(res_3<res_3(1)*1.02)
y7_3 = gsua_eval(T7_3.Estlsqc(:,1:th_3),T7_3,xdata,ydata2);
savefig('iteration4/figures/curves.fig')

bestest_3 = y7_3(1,:);
trend = [bestest_3(1),bestest_3(2:end)-bestest_3(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration4/figures/EstimatedvsReal.fig')

T7_3.Nominal = T7_3.Estlsqc(:,1);
%res : funciones de costo
th_3 = sum(res_3<res_3(1)*1.02) % threshold, cambiar el 1.01
%th = 476 %el parametro que habia dado el profe, a nosotros nos da 1????

T7_3 = gsua_ia(T7_3,T7_3.Estlsqc(:,1:th_3), false, true); 
savefig('iteration4/figures/Correlations.fig')

%Uncertainty analysis

%analisis de incertidumbre con T7_3. veamos que las curvas hagan una banda 
%al lado de los datos que estamos estimando. Si la banda es chevere,
%terminamos, si no, fijar algun parametro. 

Ua = gsua_ua(M, T7_3, 'parallel', false, 'ynom',ydata2);
savefig('iteration4/figures/Montecarlo.fig')