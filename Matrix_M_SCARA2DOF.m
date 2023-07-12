function M = Matrix_M_2DOF(m1,m2,l1,l2,q1,q2)
 
m11  = (m2*(2*l1^2 + 4*cos(q2)*l1*l2 + 2*l2^2))/2 + l1^2*m1;
m12  = l2*m2*(l2 + l1*cos(q2));

m21  = l2*m2*(l2 + l1*cos(q2));
m22  = l2^2*m2;


M = [m11 m12;
     m21 m22];
end
