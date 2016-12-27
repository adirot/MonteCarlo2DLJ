function B2 = calc_B2_hcr(T,m)
        beta = 1/T;
fun = @(z,t,n) exp(-z*t)./t.^n;
my_expint = @(n,z) integral(@(t)fun(z,t,n),1,Inf);
switch m
     case 3
        B2 = pi*(-2 + (1/3)*(6 - 4 * my_expint(5/3,beta/8) + gamma(-2/3)*(-beta)^(2/3) ));
end
        
end
