function [outputArg1,outputArg2] = iteration(inputArg1,inputArg2)
%%
%3.5. Parameter estimation

%(takes a lot of time)

opt = optimoptions('lsqcurvefit','UseParallel',true,'Display','iter');
[T7,res] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); %mover N ( 100 buena, 50 aceptable, 10 :( )
save('iteration1/values/Results7.mat','T7','res','xdata','ydata2')

%%
%3.6 Identifiability analysis 
%Cogemos el 30 % para solo tomar una de las familias
th = sum(res<res(1)*1.3)
y7 = gsua_eval(T7.Estlsqc(:,1:th),T7,xdata,ydata2);
savefig('iteration1/figures/curves.fig')

%%
bestest = y7(1,:);
trend = [bestest(1),bestest(2:end)-bestest(1:end-1)];
plot(trend,'b')
hold on
plot(diff(ydata2),'r')
title('Estimated vs Real Weekly Infections')
xlabel('Weeks')
ylabel('Cases')
legend({'Estimated','Real'})
savefig('iteration1/figures/EstimatedvsReal.fig')

%%
T7.Nominal = T7.Estlsqc(:,1);
%res : funciones de costo
th = sum(res<res(1)*1.3) % threshold, cambiar el 1.01
%th = 476 %el parametro que habia dado el profe, a nosotros nos da 1????

T7_ia = gsua_ia(T7,T7.Estlsqc(:,1:th), false, true); 
savefig('iteration1/figures/CorrelationsDiag.fig')
end

