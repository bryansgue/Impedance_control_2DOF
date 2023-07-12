function [x_sig_RK4] = xp_RK4(x,T,ts,ff,m1,m2,l1,l2,Fa)

k1= funcion_RK4(x,T,ff,m1,m2,l1,l2,Fa);
k2= funcion_RK4(x+ts/2*k1,T,ff,m1,m2,l1,l2,Fa);
k3= funcion_RK4(x+ts/2*k2,T,ff,m1,m2,l1,l2,Fa);
k4= funcion_RK4(x+ts*k3,T,ff,m1,m2,l1,l2,Fa);
x_sig_RK4= x+ts/6*(k1 +2*k2 +2*k3 +k4);
%qp_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);


%     k1 = f(st, con);   % new 
%     k2 = f(st + ts/2*k1, con); % new
%     k3 = f(st + ts/2*k2, con); % new
%     k4 = f(st + ts*k3, con); % new
%     st_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);

end