function [st_next_RK4]=x_p_RK4_2DOF(velReal,ts)

    k1 = velReal;
    k2 = velReal + ts/2*k1;% new
    k3 = velReal + ts/2*k2; % new
    k4 = velReal + ts*k3; % new
    st_next_RK4=(ts/6*(k1 +2*k2 +2*k3 +k4));
    
%     k1 = f(st, con);   % new 
%     k2 = f(st + ts/2*k1, con); % new
%     k3 = f(st + ts/2*k2, con); % new
%     k4 = f(st + ts*k3, con); % new
%     st_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);
    
end