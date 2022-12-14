function expr = x7m(t,in2,in3)
%X7M
%    EXPR = X7M(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    03-Nov-2022 16:21:25

He = in2(3,:);
Hi = in2(2,:);
Hr = in2(5,:);
Hs = in2(6,:);
Me = in2(4,:);
Mi = in2(8,:);
Ms = in2(7,:);
param64 = in3(:,1);
param65 = in3(:,2);
param66 = in3(:,3);
t2 = He.*param66;
t3 = Hi.*param64;
t4 = Me+Mi+Ms;
t5 = He.*Me.*2.0;
t6 = He.*Mi.*2.0;
t7 = He.*Ms.*2.0;
t8 = He+Hi+Hr+Hs;
t10 = Hs.*Mi.*7.311e+3;
t9 = 1.0./t4;
t11 = 1.0./t8;
t12 = -t10;
mt1 = [t2;Hi.*(-4.0e-4)+t2-t3;t9.*(t5+t6+t7+t12+Me.*t2.*5.0e+3+Mi.*t2.*5.0e+3+Ms.*t2.*5.0e+3).*(-2.0e-4);t11.*(He.*Me.*9.106e+3+Hi.*Me.*9.106e+3+Hr.*Me.*9.106e+3+Hs.*Me.*9.106e+3-Hi.*Ms.*9.93e+2).*(-1.0e-4);Hr.*(-4.0e-4)+t3;(t9.*(t5+t6+t7+t12+Hi.*Me.*2.0+Hi.*Mi.*2.0+Hr.*Me.*2.0+Hr.*Mi.*2.0+Hi.*Ms.*2.0+Hr.*Ms.*2.0))./5.0e+3];
mt2 = [(t11.*(He.*param65.*1.0e+4+Hi.*param65.*1.0e+4+Hr.*param65.*1.0e+4+Hs.*param65.*1.0e+4-He.*Ms.*1.6e+3-Hi.*Ms.*2.593e+3-Hr.*Ms.*1.6e+3-Hs.*Ms.*1.6e+3))./1.0e+4;Me.*7.506e-1-Mi.*(4.0./2.5e+1)];
expr = [mt1;mt2];
