function u = exp_bump(x)

u = zeros(size(x));
temp = 1-x.^2;
u(temp>0) = exp(1-1./temp(temp>0));

end