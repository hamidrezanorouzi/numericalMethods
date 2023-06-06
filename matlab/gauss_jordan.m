function [x] = gauss_jordan(A, b)
% A: matrix of coefficients
% b: Vertical vector of constants
% x: solution vector

A=input('enter matrix of coefficients :');
b=input('enter vector of Vertical vector :');
n = size(A,1);
Ab = [A b];

for i= 1 : n-1
    for j = i+1 : n
        Ab(j,:) = Ab(j,:) - Ab(i,:)*(Ab(j,i)/Ab(i,i));
    end
end

for i= n :-1: 2
    for j= i-1 :-1: 1
        Ab(j,:) = Ab(j,:) - Ab(i,:)*(Ab(j,i)/Ab(i,i));
    end
end

x=Ab(:,n+1)

end






