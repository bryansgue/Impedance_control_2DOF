function [s_p] = funcion_s_RK4(s,T,f,ff,m1,m2,l1,l2,q,qp)

x = s(1:2)
xp = s(3:4);

Mq = Matrix_M_SCARA2DOF(m1,m2,l1,l2,q(1),q(2));
Cq = Matrix_C_SCARA2DOF(m1,m2,l1,l2,q(1),q(2),qp(1),qp(2));
Gq = Matrix_G_SCARA2DOF(m1,m2,l1,l2,q(1),q(2));

Ja = Jaco_Scara2DOF(l1,l2,q(1),q(2));

Ja_p= JacoP_Scara2DOF(l1,l2,q(1),q(2),qp(1),qp(2));

% Conversion del modelo din�mico a coordenadas cartesianas
Mcart= inv(Ja)'*Mq*inv(Ja);
Ccart= inv(Ja)'*Cq*inv(Ja) - inv(Ja)'*Mq*inv(Ja)*Ja_p*inv(Ja);
Gcart= inv(Ja)'*Gq;

s_p = [xp;inv(Mcart)*(T-Ccart*xp-Gcart-f-ff)];


end