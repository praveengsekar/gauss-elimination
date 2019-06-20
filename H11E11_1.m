
























%code to implement gauss elimination with and without pivoting functions 

clear all
clc

for n = 2:100
    A = rand(n,n); 
    %x = rand(n,1);
    b = rand(n,1); 
    
    x = gaussel(A,b);
    y= partialpivotinggaussel(A,b);
    z = A\b;
    
    deltax(n-1) = max ( abs (x-z) );
    deltay(n-1) = max ( abs (y-z) );
end


p1 = semilogy([2:100],deltax,'b*');
hold on
p2 = semilogy([2:100],deltay,'ro');
title('Numerical Error for Gaussian eleimination with and without pivoting')
xlabel('Various sets of Equations') % x-axis label
ylabel('Numerical Error') % y-axis label
legend([p1 p2],'Without pivoting','With pivoting')

clear all

n = 50;
A = rand(n,n);
b = rand(n,1);

x = gaussel(A,b);
y = partialpivotinggaussel(A,b);
z = A\b;
deltax(1) = max ( abs (x-z) );
deltay(1) = max ( abs (y-z) );

for k = 1:9
    x = iterativegaussel(A,b,k);
    y= iterativepartialpivotinggaussel(A,b,k);
    z = A\b;
    
    deltax(k+1) = max ( abs (x-z) );
    deltay(k+1) = max ( abs (y-z) );
end

figure
p1 = semilogy(deltax,'b');
hold on
p2 = semilogy(deltay,'r');
title('Numerical Errors after iterative improvement method')
xlabel('Number of iterations') % x-axis label
ylabel('Numerical Error') % y-axis label
legend([p1 p2],'Without pivoting','With pivoting')


for k = 1:9
    if ( deltax(k) < deltax(k+1) )
        break;
    end
end
fprintf (' \n Meaningful iterations without pivoting = %d', k-1);
fprintf (' \n The error after improvement is = %d', deltax(k));


for k = 1:9
    if ( deltay(k) < deltay(k+1) )
        break;
    end
end
fprintf (' \n Meaningful iterations with pivoting = %d', k-1);
fprintf (' \n The error after improvement is = %d', deltax(k));

