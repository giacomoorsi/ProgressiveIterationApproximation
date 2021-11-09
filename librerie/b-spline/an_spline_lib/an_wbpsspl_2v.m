function [bsuv,bxsuv,bysuv,bx2suv,by2suv,bxysuv,X,Y]=an_wbpsspl_2v(gu,gv,u,v,W,x,y)
% Valutazione delle funzioni RB-spline base bivariate e delle derivate
% parziali prime e seconde in un punto (x,y);
% se x ed y sono vettori torna delle matrice 4D (per es. bsuv(l,k,:,:) si riferisce
% al punto griglia di indici (l,k), mentre bsuv(:,:,i,j) si riferisce alla
% RB-spline di indici (i,j)) e contiene i valori delle funzioni base 
% calcolate nei punti;
% bsuv contiene i valori, bxsuv contiene i valori delle derivate parziali
% prima in x, bysuv contiene i valori delle derivate parziali prima in y,
% bx2suv contiene i valori delle derivate parziali seconda in x, bysuv
% contiene i valori delle derivate parziali seconda in y e infine bxysuv
% contiene i valori delle derivate miste in x ed y
% function [bsuv,bxsuv,bysuv,bx2suv,by2suv,bxysuv,X,Y]=an_wbpsspl_2v(gu,gv,u,v,W,x,y)% gu --> grado in u della rational spline
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
% bx2suv<--matrice delle derivate parziali seconde in x
% by2suv<--matrice delle derivate parziali seconde in y
% bxysuv<--matrice delle derivate parziali miste in x ed y
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
bx2suv=zeros(nx,ny,nph,mpk);
by2suv=zeros(nx,ny,nph,mpk);
bxysuv=zeros(nx,ny,nph,mpk);

%Valutazione RB-spline bivariate e loro derivate prime parziali
 [zx,zx1,zx2]=an_bpsspl(gu,u,x);
 [zy,zy1,zy2]=an_bpsspl(gv,v,y);
 sbs=zeros(nx,ny);
 sbsx=zeros(nx,ny);
 sbsy=zeros(nx,ny);
 sbsx2=zeros(nx,ny);
 sbsy2=zeros(nx,ny);
 sbsxy=zeros(nx,ny);
 for i=1:nph
    for j=1:mpk
       bsuv(:,:,i,j)=W(i,j).*zx(:,i)*zy(:,j)';
       sbs=sbs+bsuv(:,:,i,j);
       bxsuv(:,:,i,j)=W(i,j).*zx1(:,i)*zy(:,j)';
       sbsx=sbsx+bxsuv(:,:,i,j);
       bysuv(:,:,i,j)=W(i,j).*zx(:,i)*zy1(:,j)';
       sbsy=sbsy+bysuv(:,:,i,j);
       bx2suv=W(i,j).*zx2(:,i)*zy(:,j)';
       sbsx2=sbsx2+bx2suv;
       by2suv=W(i,j).*zx(:,i)*zy2(:,j)';
       sbsy2=sbsy2+by2suv;
       bxysuv=W(i,j).*zx1(:,i)*zy1(:,j)';
       sbsxy=sbsxy+bxysuv;
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
       bx2suv(:,:,i,j)=W(i,j).*zx2(:,i)*zy(:,j)';
       bx2suv(:,:,i,j)=(bx2suv(:,:,i,j) - 2.*bxsuv(:,:,i,j).*sbsx - bsuv(:,:,i,j).*sbsx2)./sbs;
       by2suv(:,:,i,j)=W(i,j).*zx(:,i)*zy2(:,j)';
       by2suv(:,:,i,j)=(by2suv(:,:,i,j) - 2.*bysuv(:,:,i,j).*sbsy - bsuv(:,:,i,j).*sbsy2)./sbs;
       bxysuv(:,:,i,j)=W(i,j).*zx1(:,i)*zy1(:,j)';
       bxysuv(:,:,i,j)=(bxysuv(:,:,i,j) - (bysuv(:,:,i,j).*sbsy + bxsuv(:,:,i,j).*sbsx) - bsuv(:,:,i,j).*sbsxy)./sbs;
    end
 end