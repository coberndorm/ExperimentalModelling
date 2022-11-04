%Tercera Iteracion
%5.1 Aumentamos los rangos de thetaH
T.Range('lambda',1) = 1100;
T.Range('theta_h',1) = 0.2;
T.Range('Ms0',2) = 1.6e6;
T.Range('theta_m',2) = 1.2;

%%
%5.2 estimation
opt = optimoptions('lsqcurvefit','UseParallel',false,'Display','iter');
[T7_2,res_2] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); 
save('iteration3/values/Results7.mat','T7_2','res_2','xdata','ydata2')

%%
%5.2 Identifiability analysis 
th_2 = sum(res_2<res_2(1)*1.01)
y7_2 = gsua_eval(T7_2.Estlsqc(:,1:th_2),T7_2,xdata,ydata2);
savefig('iteration3/figures/curves.fig')

%%
bestest_2 = y7_2(1,:);
trend = [bestest_2(1),bestest_2(2:end)-bestest_2(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration3/figures1/EstimatedvsReal.fig')

%%
T7_2.Nominal = T7_2.Estlsqc(:,1);

T7_2 = gsua_ia(T7_2,T7_2.Estlsqc(:,1:th_2), false, true); 
savefig('iteration3/figures/Correlations.fig')
save('iteration3/values/T.mat','T')

%%
%4.4 Uncertainty analysis

%analisis de incertidumbre con T7_2. veamos que las curvas hagan una banda 
%al lado de los datos que estamos estimando. Si la banda es chevere,
%terminamos, si no, fijar algun parametro. 

Ua = gsua_ua(M, T7_2, 'parallel', false, 'ynom',ydata2);
savefig('iteration3/figures/Montecarlo.fig')
