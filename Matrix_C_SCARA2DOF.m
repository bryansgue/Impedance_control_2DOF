function M = Matrix_C_2DOF(m1,m2,l1,l2,q1,q2,q1_p,q2_p)
 
c11  = -l1*l2*m2*q2_p*sin(q2);
c12  = -l1*l2*m2*sin(q2)*(q1_p + q2_p);

c21  = l1*l2*m2*q1_p*sin(q2);
c22  = 0;


M = [c11 c12;
     c21 c22];
end