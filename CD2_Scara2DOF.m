function CD2 = CD2_Scara2DOF(l1,l2,q1,q2)

%q1 = -q1;
hx =  l2*cos(q1 + q2) + l1*cos(q1);
hy =  l2*sin(q1 + q2) + l1*sin(q1);



CD2 = [hx;hy];
return