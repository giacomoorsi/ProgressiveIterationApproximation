function y=zfunpp12(x)
%*************************************************
% derivata seconda della funzione test 
% y=(x-1)*exp(-n*x)+x^n per n=1,5,10 def. in [0,1]
%*************************************************
  n=5;
  y=n.*x.^(n - 2).*(n - 1) - 2.*n.*exp(-n.*x) + n.^2.*exp(-n.*x).*(x - 1);
end
