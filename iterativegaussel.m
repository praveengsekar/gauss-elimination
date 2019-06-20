%
%  file gaussel.m
%
%  x=gaussel(A,b) returns the solution of the set of
%  linear algebraic equations A*x=b using Gauss elimination
%  without pivoting. b and x can have several columns
%
%
function [x] = iterativegaussel(A,b,iter)

x = gaussel(A,b);

for i = 1:iter
    btilde = A*x;
    deltab = btilde - b;
    deltax = gaussel(A,deltab);
    x = x - deltax;
end


