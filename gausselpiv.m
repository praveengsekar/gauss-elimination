function [x] = gausselpiv(A,b)
%
%
% gauss elimination with partial pivoting
% input: A, matrix of system of size nxn
% input: b, matrix of free terms of size nxm
% output: x, solution of size nxm
% 
%

% first check the dimensions
[m,n] = size(A);
x=zeros(size(b)); % predefinition of x
if m~=n | n~=size(b,1) | size(x,2)~=size(b,2)
     error('not a square matrix problem')
end

B = [A b];
N = size(B,2);
% permutation vector
P = [1:n]';

for k=1:n-1    % loop over columns where the zeros will appear
    max = 0;
    index = 0;
    for j = k:n
        if abs(B(P(j),k))>max
            max = abs(B(P(j),k));
            index = j;
        end
    end
    if P(k) ~= P(index)
        temp = P(k);
        P(k) = P(index);
        P(index) = temp;
    end
    % instead of getting row k like before, now we take the row that 
    % P(k) corresponds at k
    fac=1/B(P(k),k);
    % use P(i) instead of i
    for i=k+1:n   % loop over rows where subtractions take place
        fac1=fac*B(P(i),k); % factor
        B(P(i),k)=0; % new zero by construction
        B(P(i),k+1:N)=B(P(i),k+1:N)-B(P(k),k+1:N)*fac1; % subtraction
    end
end

% Solution by backsubstitution :
% use P(k) instead of k at B matrix
for k=n:-1:1
  x(k,:)=B(P(k),n+1:N);
  for j=k+1:n
    x(k,:)=x(k,:)-B(P(k),j)*x(j,:);
  end
  x(k,:)=x(k,:)/B(P(k),k);
end
