clear all
clc

n = 50;
m = 1;
A = rand(n,n);
b = rand(n,m);

[x,c_el,c_sub] = gausselcounter(A,b);
f_el = 0;
for i = 1:n
    f_el = f_el + ( (i-1) * ( i+m) );
end

f_sub = m * n * (n+1) * (1/2);

fprintf('n = %d \n m = %d', n,m);
fprintf('\nThe number of multiplication calculated from the gaussel function');
fprintf('\n The number of mul/div in elimination is %d', c_el);
fprintf('\n The number of mul/div in substitution is %d', c_sub);

fprintf('\n\nThe number of multiplication calculated using formula');
fprintf('\n The number of mul/div in elimination is %d', f_el);
fprintf('\n The number of mul/div in substitution is %d', f_sub);

n = 50;
m = 2;
A = rand(n,n);
b = rand(n,m);

[x,c_el,c_sub] = gausselcounter(A,b);
f_el = 0;
for i = 1:n
    f_el = f_el + ( (i-1) * ( i+m) );
end

f_sub = m * n * (n+1) * (1/2);

fprintf('\n\nn = %d \n m = %d', n,m);
fprintf('\nThe number of multiplication calculated from the gaussel function');
fprintf('\n The number of mul/div in elimination is %d', c_el);
fprintf('\n The number of mul/div in substitution is %d', c_sub);

fprintf('\n\nThe number of multiplication calculated using formula');
fprintf('\n The number of mul/div in elimination is %d', f_el);
fprintf('\n The number of mul/div in substitution is %d', f_sub);
