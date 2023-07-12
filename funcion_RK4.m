function [x_p] = funcion_RK4(x,T,ff,m1,m2,l1,l2,Fa)

q = x(1:2);
qp = x(3:4);

Ja = Jaco_Scara2DOF(l1,l2,q(1),q(2));

Mq = Matrix_M_SCARA2DOF(m1,m2,l1,l2,q(1),q(2));
Cq = Matrix_C_SCARA2DOF(m1,m2,l1,l2,q(1),q(2),qp(1),qp(2));
Gq = Matrix_G_SCARA2DOF(m1,m2,l1,l2,q(1),q(2));

x_p = [qp;inv(Mq)*(T+Ja'*Fa-Cq*qp-Gq-ff)];


end