function [cpr,tr]=an_knotremMQ(g,cp,t,tc)
% Rimuove i nodi tc dalla partizione estesa t, e determina i coefficienti
% nella nuova rappresentazione, mediante soluzione di un sistema 
% sovradeterminato
% function [cpr,tr]=an_knotrem(g,cp,t,tc)
% g   --> grado della spline
% cp  --> vettore dei coefficienti
% t   --> vettore dei nodi partizione estesa
% tc  --> vettore dei nodi da rimuovere
% cpr <-- vettore dei nuovi coefficienti
% tr  <-- vettore nuovo dei nodi partizione estesa
nt=length(t);
ntc=length(tc);
ncp=length(cp);
n=nt-g-1;
M=eye(ncp,ncp);

   for l=1:ntc
      nt=nt-1;
      ncp=ncp-1;
      A=zeros(ncp+1,ncp);
%determina intervallo nodale contenente tc(l)
%cerca intervallo (t(j),t(j+1)]
%while (tc(l)>=t(j+1)) e' ugualmente corretto; cerco l'intervallo
%destro contenente tc(l), anziche' il sinistro
      j=0;
      while(tc(l)>t(j+1))
        j=j+1;
      end
%costruisce matrice rettangolare (ncp+1)x(ncp)
      for i=1:j-g
         A(i,i)=1.0;
      end
      for i=j-g+1:j
         lambda=(tc(l)-t(i))/(t(i+g+1)-t(i));
         A(i,i-1)=1-lambda;
         A(i,i)=lambda;
      end
      for i=j+1:ncp+1
         A(i,i-1)=1.0;
      end
      M=M*A;
%       disp(A);
%       disp(rank(A));
%rimozione del nodo tc(l) dalla partizione
      for k=j+1:nt
        t(k)=t(k+1);
      end
      clear A;
  end
% disp('matrice rettangolare del sistema sovradeterminato');
% disp(M);
% risolve sistema sovradeterminato
cpr=(M\cp')';
tr=t(1:nt);