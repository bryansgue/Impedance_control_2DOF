function [st_next_RK4]=RK4_2DOF(velReal,Mcart, Ccart, Gcart,ft, ff, f,ts)

    k1 = DinInversa2DOF(velReal,Mcart, Ccart, Gcart,ft, ff, f);
    k2 = DinInversa2DOF(velReal + ts/2*k1, Mcart, Ccart, Gcart,ft, ff, f);% new
    k3 = DinInversa2DOF(velReal + ts/2*k2, Mcart, Ccart, Gcart,ft, ff, f); % new
    k4 = DinInversa2DOF(velReal + ts*k3, Mcart, Ccart, Gcart,ft, ff, f); % new
    st_next_RK4=(ts/6*(k1 +2*k2 +2*k3 +k4));
    
%     k1 = f(st, con);   % new 
%     k2 = f(st + ts/2*k1, con); % new
%     k3 = f(st + ts/2*k2, con); % new
%     k4 = f(st + ts*k3, con); % new
%     st_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);
    
end