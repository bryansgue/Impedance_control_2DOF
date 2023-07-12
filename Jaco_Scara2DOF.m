function J = Jaco_Scara2DOF(l1,l2,q1,q2)


J11 = - l2*sin(q1 + q2) - l1*sin(q1);
J12 = -l2*sin(q1 + q2);

J21 = l2*cos(q1 + q2) + l1*cos(q1);
J22 = l2*cos(q1 + q2);


J = [J11 J12;
     J21 J22];

end