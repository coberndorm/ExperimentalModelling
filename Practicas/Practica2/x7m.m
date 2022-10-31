function expr = x7m(t,in2,in3)
%X7M
%    EXPR = X7M(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    30-Oct-2022 23:00:24

He = in2(3,:);
Hi = in2(2,:);
Hr = in2(4,:);
Hs = in2(6,:);
Me = in2(5,:);
Mi = in2(8,:);
Ms = in2(7,:);
param78 = in3(:,1);
param79 = in3(:,2);
param80 = in3(:,3);
param81 = in3(:,4);
param82 = in3(:,5);
param83 = in3(:,6);
param84 = in3(:,7);
t2 = He.*param83;
t3 = Hi.*param80;
t4 = He.*Me;
t5 = He.*Mi;
t6 = He.*Ms;
t7 = Me+Mi+Ms;
t8 = Hi.*Ms.*param79;
t9 = He+Hi+Hr+Hs;
t11 = Hs.*Mi.*param78.*2.5e+3;
t10 = 1.0./t7;
t12 = 1.0./t9;
t13 = -t11;
expr = [t2;Hi.*(-4.0e-4)+t2-t3;t10.*(t4+t5+t6+t13+Me.*t2.*2.5e+3+Mi.*t2.*2.5e+3+Ms.*t2.*2.5e+3).*(-4.0e-4);Hr.*(-4.0e-4)+t3;-t12.*(-t8+param82.*t4+param84.*t4+Hi.*Me.*param82+Hi.*Me.*param84+Hr.*Me.*param82+Hr.*Me.*param84+Hs.*Me.*param82+Hs.*Me.*param84);(t10.*(t4+t5+t6+t13+Hi.*Me+Hi.*Mi+Hr.*Me+Hr.*Mi+Hi.*Ms+Hr.*Ms))./2.5e+3;-t12.*(t8-He.*param81-Hi.*param81-Hr.*param81-Hs.*param81+param82.*t6+Hi.*Ms.*param82+Hr.*Ms.*param82+Hs.*Ms.*param82);Me.*param84-Mi.*param82];