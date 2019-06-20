%
%  file gaussel.m
%
%  x=gaussel(A,b) returns the solution of the set of
%  linear algebraic equations A*x=b using Gauss elimination
%  without pivoting. b and x can have several columns
%
%
function [x] = pgaussel(A,b)
[m,n]=size(A);
if m~=n | n~=size(b,1), error('not a square matrix problem'); end;

B=[A b];
N=size(B,2);
p = 1:m;
for i = 1:m-1
        maximum = max( B(p(i:m) , i) );
        index = i-1+find( B(p(i:m),i) == maximum) ;
        temp = p(index);
        p(index) = p(i);
        p(i) = temp;
end

% bring the matrix into triangular form (Gauss elimination):

for k=1:n-1,    % loop over columns where the zeros will appear
  fac=1/B(p(k),k);
  for i=k+1:n   % loop over rows where subtractions take place
    fac1=fac*B(p(i),k); % factor
    B(p(i),k)=0; % new zero by construction
    B(p(i),k+1:N)=B(p(i),k+1:N)-B(p(k),k+1:N)*fac1; % subtraction
  end
end

% Solution by backsubstitution :
x=zeros(size(b)); % predefinition of x
for k=n:-1:1
  x(k,:)=B(p(k),n+1:N);
  for j=k+1:n
    x(k,:)=x(k,:)-B(p(k),j)*x(j,:);
  end
  x(k,:)=x(k,:)/B(p(k),k);
end