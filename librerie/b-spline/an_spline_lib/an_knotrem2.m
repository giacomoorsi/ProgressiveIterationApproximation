function [cpr,tr]=an_knotrem2(g,cp,t,tc)
% Rimuove i nodi tc dalla partizione estesa t, e determina i coefficienti
% nella nuova rappresentazione
% function [cpr,tr]=an_knotrem2(g,cp,t,tc)
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

   for l=1:ntc
      nt=nt-1;
      ncp=ncp-1;

      j=0;
%determina intervallo nodale contenente tc(l)
%cerca intervallo (t(j),t(j+1)]
%while (tc(l)>=t(j+1)) e' ugualmente corretto; cerco l'intervallo
%destro contenente tc(l), anziche' il sinistro
      while(tc(l)>t(j+1))
        j=j+1;
      end
%aggiornamneto dei CP; da 1:j-g restano gli stessi
      k=j-g+1;
      while (k<=j)
        lambda=(tc(l)-t(k))/(t(k+g+1)-t(k));
        cp(k)=(cp(k)-(1-lambda)*cp(k-1))/lambda;
        k=k+1;
      end
%controllo su rimozione; deve essere cp(j)==cp(j+1)
% if(abs(cp(j)-cp(j+1)) < 1.0e-5)
%     disp('rimuovibile')
% else
%     disp('non rimuovibile in modo esatto')
% end
      for k=j+1:ncp
        cp(k)=cp(k+1);
      end
%rimozione del nodo nella partizione */
      for k=j+1:nt
        t(k)=t(k+1);
      end
    end
tr=t(1:nt);
cpr=cp(1:ncp);

