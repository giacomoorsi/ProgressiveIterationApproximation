function bs=an_bspl_sym(g,t,el,id)
%Determina l'espressione simbolica delle B-spline di grado g (definite 
%su una partizione nodale) non nulle nell'intervallo [t(el), t(el+1)];
%deriva dal codice gc_bspl.m di cui segue la stessa implementazione
% function bs=an_bspl_sym(g,t,el,id)
% g  --> grado della spline
% t  --> vettore dei nodi partizione estesa
% el --> indice dell'intervallo non degenere, cioe' [t(el),t(el+1)) con
%        t(el) < t(el+1)
% id --> id<0 non stampa nulla
%        id=0 stampa tutte le B-spline non nulle sull'intervallo [t(el),t(el+1))
%        id=1 stampa solo la prima non nulla sull'intervallo, 
%        id=2 la seconda non nulla sull'intervallo,
%        ecc.
% bs <-- vettore simbolico delle funzioni B-spline nell'intervallo
syms x
syms temp d1 d2 gcbeta

if (nargin<4)
    id=-1;
end

dim=length(t)-g-1;
bs=zeros(1,dim);
bs=sym(bs);
if (t(el) < t(el+1))
    bs(el)=1;
    k=el;
    for i=1:g
       temp=0;
       for j=el:k
         d1=x-t(j);
         d2=t(i+j)-x;
         gcbeta=bs(j)/(d1+d2);
         bs(j-1)=temp+d2*gcbeta;
         temp=d1*gcbeta;
       end
       bs(k)=temp;
       el=el-1;
    end

    if (id>=0)
      if (id==0)
        for j=k-g:k
        disp('***************************************************');
        pretty(simplify(eval(bs(j))));
        end
%         disp('********************************');
%         disp('verifica simbolica di somma ad 1');
%         disp(simplify(sum(bs)));
      else
        if (1<= id & id<=g+1) 
          j=k-g+id-1;
          pretty(simplify(eval(bs(j))));
        else
          disp('indice non lecito');
        end
      end
    end
else
    disp('intervallo nodale degenere');
end
