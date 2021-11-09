function [a, b] = test_fun(i)
 switch (i)
     case 'fun1'
         a = -pi;
         b = pi;
     case 'fun2' 
         a = -2;
         b = 2;
     case 'fun3'
         a = -2;
         b = 1;
     case 'gerono';
         a = 0;
         b = 2*pi;
         
     otherwise 
         a = 0;
         b = 0;
end