function I = getLoadCurrent(t,t0,tf);

I_final = 1;
tau = 100;
I = I_final*(1- exp(-t/tau));