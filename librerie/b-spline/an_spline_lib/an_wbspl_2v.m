function [bsuv,X,Y]=an_wbspl_2v(gu,gv,u,v,W,x,y)
% Valutazione delle funzioni RB-spline base bivariate in un punto (x,y);
% se x ed y sono vettori torna la matrice bsuv 4D (dove bsuv(l,k,:,:) si riferisce
% al punto griglia di indici (l,k), mentre bsuv(:,:,i,j) si riferisce alla
% RB-spline di indici (i,j)) e contiene i valori delle funzioni base 
% calcolate nei punti;
% function [bsuv,X,Y]=an_wbspl_2v(gu,gv,u,v,W,x,y)
% gu --> grado in u della rational spline
% gv --> grado in v della rational spline
% u --> vettore dei nodi partizione estesa
% v --> vettore dei nodi partizione estesa
% W --> matrice dei pesi
% x --> vettore di punti
% y --> vettore di punti
% bsuv <-- matrice 4D delle funzioni RB-spline sulla griglia di punti
%          bsuv(:,:,i,j) rappresenta i valori della RB-spline (i,j) su
%          tutti i punti della griglia di valutazione
% X <-- matrice delle ascisse
% Y <-- matrice delle ordinate
nu=length(u);
nph=nu-gu-1;
nv=length(v);
mpk=nv-gv-1;

[X,Y]=meshgrid(x,y);
X=X';
Y=Y';
nx=length(x);
ny=length(y);
bsuv=zeros(nx,ny,nph,mpk);

%Valutazione RB-spline bivariate
zx=an_bspl(gu,u,x);
zy=an_bspl(gv,v,y);
sbs=zeros(nx,ny);
for i=1:nph
   for j=1:mpk
       bsuv(:,:,i,j)=W(i,j).*zx(:,i)*zy(:,j)';
       sbs=sbs+bsuv(:,:,i,j);
   end
end
 
for i=1:nph
   for j=1:mpk
       bsuv(:,:,i,j)=W(i,j).*zx(:,i)*zy(:,j)';
       bsuv(:,:,i,j)=bsuv(:,:,i,j)./sbs;
   end
end
