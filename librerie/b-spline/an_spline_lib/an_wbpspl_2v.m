function [bsuv,bxsuv,bysuv,X,Y]=an_wbpspl_2v(gu,gv,u,v,W,x,y)
% Valutazione delle funzioni RB-spline base bivariate e delle derivate
% parziali prime in un punto (x,y);
% se x ed y sono vettori torna delle matrice 4D (per es. bsuv(l,k,:,:) si riferisce
% al punto griglia di indici (l,k), mentre bsuv(:,:,i,j) si riferisce alla
% RB-spline di indici (i,j)) e contiene i valori delle funzioni base 
% calcolate nei punti;
% bsuv contiene i valori, bxsuv contiene i valori delle derivate parziali
% prima in x, bysuv contiene i valori delle derivate parziali prima in y;
% function [bsuv,bxsuv,bysuv,X,Y]=an_wbpspl_2v(gu,gv,u,v,W,x,y)
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
% bxsuv<-- matrice delle derivate parziali prime in x
% bysuv<-- matrice delle derivate parziali prime in y
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
bxsuv=zeros(nx,ny,nph,mpk);
bysuv=zeros(nx,ny,nph,mpk);

%Valutazione RB-spline bivariate e loro derivate prime parziali
 [zx,zx1]=an_bpspl(gu,u,x);
 [zy,zy1]=an_bpspl(gv,v,y);
 sbs=zeros(nx,ny);
 sbsx=zeros(nx,ny);
 sbsy=zeros(nx,ny);
 for i=1:nph
    for j=1:mpk
       bsuv(:,:,i,j)=W(i,j).*zx(:,i)*zy(:,j)';
       sbs=sbs+bsuv(:,:,i,j);
       bxsuv(:,:,i,j)=W(i,j).*zx1(:,i)*zy(:,j)';
       sbsx=sbsx+bxsuv(:,:,i,j);
       bysuv(:,:,i,j)=W(i,j).*zx(:,i)*zy1(:,j)';
       sbsy=sbsy+bysuv(:,:,i,j);
    end
 end

 for i=1:nph
    for j=1:mpk
       bsuv(:,:,i,j)=W(i,j).*zx(:,i)*zy(:,j)';
       bsuv(:,:,i,j)=bsuv(:,:,i,j)./sbs;
       bxsuv(:,:,i,j)=W(i,j).*zx1(:,i)*zy(:,j)';
       bxsuv(:,:,i,j)=(bxsuv(:,:,i,j) - bsuv(:,:,i,j).*sbsx)./sbs;
       bysuv(:,:,i,j)=W(i,j).*zx(:,i)*zy1(:,j)';
       bysuv(:,:,i,j)=(bysuv(:,:,i,j) - bsuv(:,:,i,j).*sbsy)./sbs;
    end
 end