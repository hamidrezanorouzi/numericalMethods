%------------- Matlab ---------------
% Numerical methods course, Amirkabir University of Technology
% website: www.cemf.ir
% Root finding based on the bisection method in interval (a b)

%inputs 
a = 2;
b = 5;
tol = 1.0e-5; %relative tolerance 
nMax = 100; 

fa = fx(a);
fb = fx(b);
found = false;
c_old = a;

%main loop
for iter=1:nMax
    
    c = (a+b)/2;
    fc = fx(c);
    if( abs(fc) < tol)
        found = true;
        break;
    end
        
    if( fc*fb < 0 )
        a = c;
        fa = fc;
    else
        b = c;
        fb = fc;
    end
    
    if( abs(c) > 1.0e-15 )
        ea = abs((c-c_old)/c);
    else
        ea = abs(c);
    end
    
    if( ea < tol)
       found = true;
       break;
    end
    c_old = c;
end

%display results
if( found ) 
   fprintf( 'A root has been found : %f\n', c);
   fprintf( 'Function value at root : %f\n', fc);
   fprintf( 'Number of iterations  : %d\n', iter);
else
   fprintf( 'No root was found in %d iterations\n', nMax);
end

