function CD1 = CD1_Scara2DOF(l1,l2,q1,q2)

%q1 = -q1;
hx = l1*cos(q1);
hy =  l1*sin(q1);



CD1 = [hx;hy];
return