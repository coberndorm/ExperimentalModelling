%%
%1. Define equations for the model and load timeseries data

syms  Ms(t) Me(t) Mi(t) Hs(t) He(t) Hi(t) Hr(t) Hit(t)

syms lambda beta_m mu_m theta_m mu_h beta_h theta_h...
    gamma_h

%pulse(t)=1;
H=Hs+He+Hi+Hr;
M=Ms+Me+Mi;
ode1 = diff(Ms) == lambda - beta_m*Hi*Ms/H - (mu_m)*Ms;
ode2 = diff(Me) == beta_m*Hi*Ms/H - (theta_m+mu_m)*Me;
ode3 = diff(Mi) == theta_m*Me - mu_m*Mi;
ode4 = diff(Hs) == -beta_h*Mi*Hs/M + (He+Hi+Hr)*mu_h;
ode5 = diff(He) == beta_h*Mi*Hs/M - (theta_h+mu_h)*He;
ode6 = diff(Hi) == theta_h*He - (gamma_h+mu_h)*Hi;
ode7 = diff(Hr) == gamma_h*Hi - mu_h*Hr;
ode8 = diff(Hit) == theta_h*He;
odes=[ode1; ode2 ;ode3; ode4; ode5 ;ode6; ode7; ode8];
vars=[Hit Hi Me Hr Hs He Ms Mi];
opts = odeset('NonNegative',1:8);

%2. Data preparation

load Range7.mat
load('DataBello_full.mat');
ydata=DataBello.cases(153:233)';
xdata=linspace(0,length(ydata)-1,length(ydata));
ydata2=ydata;
for i=1:length(ydata)
ydata2(i)=sum(ydata(1:i));
end

%%
%3.Analisis de sensibilidad

%3.1. Data prep

[T,~]=gsua_dpmat(odes,vars,[0 80],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output = 1;
M = gsua_dmatrix(T,100);

ynom = gsua_eval(T.Nominal,T);

plot(cumsum(ydata),'b')
title('Original Acumulated Human Infections (Hit)')
xlabel('Weeks')
ylabel('Cases')
savefig('iteration1/figures/NominalValues.fig')

%%


%3.2. Sensibility analysis

Tsa = gsua_sa(M,T,'parallel', false, 'SensMethod', 'Xiao', 'ynom', ynom);

gsua_plot('Bar',Tsa,Tsa.STi)
savefig('iteration1/figures/SensibilityAnalisis.fig')

%3.3. Confiabilidad

c = sum(Tsa.Si)/sum(abs(Tsa.Si))

%%
%3.4 Fixing parameters

%According to the sensitivity analysis, we will fix the 4 parameters with
%the lower sensitivity indexes (with the nominal values).

%These are He0
%vars=[Hit Hi Me Hr Hs He Ms Mi];
vars=[Hit Hi He Me Hr Hs Ms Mi];

HeO = T.Nominal('He0');

RangeTemp = Range; 
Range(3,:) = HeO; 
Range(4:6,:) = [RangeTemp(3,:); RangeTemp(4,:); RangeTemp(5,:)];

[T,~]=gsua_dpmat(odes,vars,[0 80],'7m','output',1,'opt',opts,'Range',Range);
T.Properties.CustomProperties.output = 1;
M = gsua_dmatrix(T,100);
save('iteration1/values/T.mat','T')
%%
%3.5. Parameter estimation

%(takes a lot of time)

opt = optimoptions('lsqcurvefit','UseParallel',true,'Display','iter');
[T7,res] = gsua_pe(T,xdata,ydata2,'N',100,'opt',opt); %mover N ( 100 buena, 50 aceptable, 10 :( )
save('iteration1/values/Results7.mat','T7','res','xdata','ydata2')

%%
%3.6 Identifiability analysis 
%Cogemos el 30 % para solo tomar una de las familias
th = sum(res<res(1)*1.01)
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
th = sum(res<res(1)*1.01) % threshold, cambiar el 1.01
%th = 476 %el parametro que habia dado el profe, a nosotros nos da 1????

T7_ia = gsua_ia(T7,T7.Estlsqc(:,1:th), false, true); 
%savefig('iteration1/figures/CorrelationsDiag.fig')



