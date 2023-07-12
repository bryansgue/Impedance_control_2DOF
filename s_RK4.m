function [x_sig_RK4] = s_RK4(s,T,f,fr,ts,m1,m2,l1,l2,q,qp)

k1= funcion_s_RK4(s,T,f,fr,m1,m2,l1,l2,q,qp);
k2= funcion_s_RK4(s+ts/2*k1,T,f,fr,m1,m2,l1,l2,q,qp);
k3= funcion_s_RK4(s+ts/2*k2,T,f,fr,m1,m2,l1,l2,q,qp);
k4= funcion_s_RK4(s+ts*k3,T,f,fr,m1,m2,l1,l2,q,qp);
x_sig_RK4= s+ts/6*(k1 +2*k2 +2*k3 +k4);
%qp_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);


%     k1 = f(st, con);   % new 
%     k2 = f(st + ts/2*k1, con); % new
%     k3 = f(st + ts/2*k2, con); % new
%     k4 = f(st + ts*k3, con); % new
%     st_next_RK4=st +ts/6*(k1 +2*k2 +2*k3 +k4);

end