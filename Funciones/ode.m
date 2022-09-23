function dxdt = ode(t, x, u, tiempo, odes)
    u_int = interp1(tiempo, u, t);
    try 
       dxdt = odes(u_int,x(1),x(2));
    catch
       dxdt = odes(x(1),x(2));
    end
end

