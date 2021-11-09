function tt = curv3_param(Q,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function tt=curv3_param(Q,param)
%Calcola la parametrizzazionedi una lista di punti 3D 
%Q      --> matrice dei punti da parametrizzare (npuntix3)
%param  --> scelta della parametrizzazione
%           0=uniforme, 1=centripeta, 2=corda
%tt     <-- parametrizzazione della lista di punti 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = Q(:,1);
Y = Q(:,2);
Z = Q(:,3);
n=length(X);
alpha=0.5*param;
tau(1)=0;
for i=2:n
  tau(i)=tau(i-1)+(sqrt((X(i)-X(i-1)).^2+(Y(i)-Y(i-1)).^2+(Z(i)-Z(i-1)).^2)).^alpha;
end
for i=1:n
  tt(i)=tau(i)/tau(n);
end