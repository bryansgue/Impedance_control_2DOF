function [h_pp] = DinInversa2DOF(velReal,Mcart, Ccart, Gcart,ft, ff, f)

h_pp = inv(Mcart)*(ft - Ccart*velReal - Gcart - ff - f);

end