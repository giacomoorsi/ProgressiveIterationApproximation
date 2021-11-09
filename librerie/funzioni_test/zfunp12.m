function y=zfunp12(x)
%*************************************************
% derivata prima della funzione test 
% y=(x-1)*exp(-n*x)+x^n per n=1,5,10 def. in [0,1]
%*************************************************
  n=5;
  y=exp(-n.*x) + n.*x.^(n - 1) - n.*exp(-n.*x).*(x - 1);
end
