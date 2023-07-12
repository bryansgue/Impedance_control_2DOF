function G = Matrix_G_2DOF(m1,m2,l1,l2,q1,q2)
 
g = 9.8;

g11  = g*(l1*m1*cos(q1) + l1*m2*cos(q1) + l2*m2*cos(q1 + q2));


g21  = g*l2*m2*cos(q1 + q2);



G = [g11;
     g21];
end