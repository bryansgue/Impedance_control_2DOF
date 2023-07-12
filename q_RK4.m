function [x_sig_RK4] = q_RK4(q,velReal,l1,l2,ts)

k1= funcion_q_RK4(q,velReal,l1,l2);
k2= funcion_q_RK4(q+ts/2*k1,velReal,l1,l2);
k3= funcion_q_RK4(q+ts/2*k2,velReal,l1,l2);
k4= funcion_q_RK4(q+ts*k3,velReal,l1,l2);
x_sig_RK4= q+ts/6*(k1 +2*k2 +2*k3 +k4);
%qp_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);


%     k1 = f(st, con);   % new 
%     k2 = f(st + ts/2*k1, con); % new
%     k3 = f(st + ts/2*k2, con); % new
%     k4 = f(st + ts*k3, con); % new
%     st_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);

end