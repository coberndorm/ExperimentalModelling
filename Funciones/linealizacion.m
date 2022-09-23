function [A,B,C,D] = linealizacion(odes,y,vars,u)
%LINEALIZACION Summary of this function goes here
%   Detailed explanation goes here

A = jacobian(odes, vars); B = jacobian(odes, u); 
C = jacobian(y,vars); 
size1 = size(B); size2 = size(C);
D = zeros(size2(1),size1(2));
end

