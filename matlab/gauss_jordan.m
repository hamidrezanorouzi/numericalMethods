function x = gauss_jordan(A,b)

% A: matrix of coefficients
% b: Vertical vector of constants
% x: solution vector

 n=size(b,1);   
 x=zeros(n,1);  
 Ab=[A b];

%First, we must check that no element on the main diameter is zero
for i=1:n-1
    for j=i+1:n
        if (Ab(i,i) == 0) 
            temp = Ab(i,:);
            Ab(i,:) = Ab(j,:);
           Ab(j,:) = temp;
        end
    end

    for i=1:n      
        Ab(i,:)=Ab(i,:)/Ab(i,i); 
        for j=1:n  
           if i ~= j  
               B= Ab(j,i);
               Ab(j,:)=Ab(j,:)-B*Ab(i,:); 
           end 
        end
    end

    x=Ab(:,n+1); 

end