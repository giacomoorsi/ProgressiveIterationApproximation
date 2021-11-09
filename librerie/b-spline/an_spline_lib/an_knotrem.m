function [cpr,tr]=an_knotrem(g,cp,t,tc)
% Rimuove i nodi tc dalla partizione estesa t, e determina i coefficienti
% nella nuova rappresentazione
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
%i CP da 1:j-g restano gli stessi
% 
%aggiornamento dei rimanenti CP
       k1=j-g+1;
       k2=j;
       while (k1<k2)
        lambda=(tc(l)-t(k1))/(t(k1+g+1)-t(k1));
        cp(k1)=(cp(k1)-(1-lambda)*cp(k1-1))/lambda;
        k1=k1+1;
        if(k1<k2)
            lambda=(tc(l)-t(k2))/(t(k2+g+1)-t(k2));
            cp(k2-1)=(cp(k2)-lambda*cp(k2+1))/(1-lambda);
            k2=k2-1;
        end
       end
%controllo su rimozione; deve essere cp(j)==cp(j+1)
%       fprintf('\ncompare %f %f \n',cp(j),cp(j+1))

%i CP da j+1:ncp vengono shiftati di indice
      for k=j:ncp
        cp(k)=cp(k+1);
      end
%rimozione del nodo nella partizione */
      for k=j+1:nt
        t(k)=t(k+1);
      end
    end
tr=t(1:nt);
cpr=cp(1:ncp);

