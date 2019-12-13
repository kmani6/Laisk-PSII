function[Fl] = LaiskFluorescence(yidcs, kidcs, k, Sol) 


yopoax = yidcs{1};
yoprao = yidcs{2};
yoprar = yidcs{3};
yrprao = yidcs{4};
yrprar = yidcs{5};

kn = kidcs(1);
kp = kidcs(2);
kr = kidcs(3); 
kq = kidcs(4);



Fl = 1/(1+k(kn)+k(kr)+k(kq))*(sum(Sol(yopoax,:)))...
     ...
     +1/(1+k(kp)+k(kn)+k(kr))*(sum(Sol(yoprao,:)))...
     ...
     +1/(1+k(kn)+k(kr))*(sum(Sol(yoprar,:)))...
     ...
     +1/(1+k(kp)+k(kn))*(sum(Sol(yrprao,:)))+...
     ...
     1/(1+k(kn))*(sum(Sol(yrprar,:)));

end
