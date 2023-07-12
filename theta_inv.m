function x_sol = theta_inv(h)
    syms x
    x_sol = vpasolve(x + (1+x)*log2(1+x)-x*log2(x) == h, x);
end

