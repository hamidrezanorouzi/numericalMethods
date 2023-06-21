A = [0 1;1 0];
b = [10; -4];

x = Gauss_Jordan_updated(A, b);

expected_x = [-4 ;10];

if isequal(x, expected_x)
    disp('The code is correct.');
else
    disp('The code is incorrect.');
end
% کد بالا صرفا برای بررسی این بود که کد درصورتی که 
% در قظر اصلی عنصری برابر صفر باشد درست عمل میکند
% : به ظور کلی میتوان ار کذ زیر نیز استفاده کرذ
% A=[3 5;2 3];
% b=[9;5];
% x = Gauss_Jordan_updated(A,b)
% expected_x = [-2 ;3];