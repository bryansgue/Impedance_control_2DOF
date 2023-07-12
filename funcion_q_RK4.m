function [q_p] = funcion_q_RK4(q,velReal,l1,l2)

Ja = Jaco_Scara2DOF(l1,l2,q(1),q(2));

q_p = inv(Ja)*velReal;

end