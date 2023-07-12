function [q_p] = h2q_2DOF(velReal,Ja)

q_p = inv(Ja)*velReal;



end