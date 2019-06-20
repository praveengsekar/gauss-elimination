%
%  file partialpivotinggaussel.m
%
%  x=partialpivotinggaussel(A,b) returns the solution of the set of
%  linear algebraic equations A*x=b using Gauss elimination
%  with partial pivoting. b and x can have several columns
%
%
function [x] = partialpivotinggaussel(A,b)
[m,n]=size(A);
d = det(A);
if m~=n | n~=size(b,1), error('not a square matrix problem'); end;
if d==0, error('not a singular matrix'); end;

B=[A b];
N=size(B,2);

% rearraging the vector p 
p = 1:m;
p = p';
for k = 1:n-1
        maximum = max( B(p(k:n) , k) ); % finding the maximum in that columns considering only the elements of that row and row below it
        index = find( B(p(1:n),k) == maximum) ; % finding the index of that maximum
        
        % swaping the positions in vector p
        if p(k) ~= p(index)
            temp = p(index); 
            p(index) = p(k);
            p(k) = temp;
        end
% bring the matrix into triangular form (Gauss elimination):

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

