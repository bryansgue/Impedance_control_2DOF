function [st_next_RK4]=q_p_RK4_2DOF(velReal,Ja,ts)

    
    k1 = h2q_2DOF(velReal,Ja);
    k2 = h2q_2DOF(velReal + ts/2*k1,Ja);
    k3 = h2q_2DOF(velReal + ts/2*k2,Ja);
    k4 = h2q_2DOF(velReal + ts*k3,Ja);
    st_next_RK4=(ts/6*(k1 +2*k2 +2*k3 +k4));
    
%     k1 = f(st, con);   % new 
%     k2 = f(st + ts/2*k1, con); % new
%     k3 = f(st + ts/2*k2, con); % new
%     k4 = f(st + ts*k3, con); % new
%     st_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);
    
end