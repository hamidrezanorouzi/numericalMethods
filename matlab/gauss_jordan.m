function x = gauss_jordan(A,b)
% A: matrix of coefficients
% b: Vertical vector of constants
% x: solution vector

A=input('enter matrix of coefficients :');
b=input('enter vector of Vertical vector :');

    n=size(b,1);   
    x=zeros(n,1);  
    Ab=[A b];
    
    for i=1:n      
        Ab(i,:)=Ab(i,:)./Ab(i,i); 
        for j=1:n  
           if i ~= j  
               Ab(j,:)=Ab(j,:)-Ab(i,:)*Ab(j,i); 
           end 
        end
    end
    x=Ab(:,n+1); 
end