function x = softnorm(xOld,q1,q2)
% Softnormalizes a signal xOld to a new signal x using q1 and q2
% If xOld is a matrix, this will be performed columnwise
    qs = quantile(xOld,[q1 q2]);    
    x = 2*(xOld-qs(1))/(qs(2)-qs(1)) - 1;
end