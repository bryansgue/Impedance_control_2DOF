function J_p = JacoP_Scara2DOF(l1,l2,q1,q2,q1_p,q2_p)

J11 = - q1_p*(l2*cos(q1 + q2) + l1*cos(q1)) - l2*q2_p*cos(q1 + q2); 
J12 = -l2*cos(q1 + q2)*(q1_p + q2_p);

J21 = - q1_p*(l2*sin(q1 + q2) + l1*sin(q1)) - l2*q2_p*sin(q1 + q2);
J22 = -l2*sin(q1 + q2)*(q1_p + q2_p);


J_p = [J11 J12;
     J21 J22];

end