function tt = curv2_param(param,X,Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function tt=curv2_param(param,X,Y)
%Calcola la parametrizzazionei di una lista di punti 2D  
%param  --> scelta della parametrizzazione
%           0=uniforme, 1=centripeta, 2=corda
%X      --> vettore delle ascisse
%Y      --> vettore delle ordinate
%tt     <-- parametrizzazione della lista di punti 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=0.5*param;
bma=max(X)-min(X);
n=length(X);
tau(1)=0;
for i=2:n
  tau(i)=tau(i-1)+(sqrt((X(i)-X(i-1)).^2+(Y(i)-Y(i-1)).^2)).^alpha;
end
for i=1:n
  tt(i)=min(X)+bma*tau(i)/tau(n);
end