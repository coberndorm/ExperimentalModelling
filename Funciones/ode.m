function dxdt = ode(t, x, u, tiempo, odes)
    u_int = interp1(tiempo, u, t);
    try 
       dxdt = odes(x(1),x(2),u_int);
    catch
       dxdt = odes(x(1),x(2));
    end
end

