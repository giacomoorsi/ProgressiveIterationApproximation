% Script main_spline_curve_ki_de
% Questo script permette di modellare una curva spline 2D
% si danno i punti interattivamente con button1; dato l'ultimo punto
% si prema button3 per fine inserimento e curva aperta, oppure button2
% per fine inserimento e curva chiusa.
% Per uscire dalla fase di modellazione dare button1 sull'area rettangolare
% grigia esterna all'area di disegno.
clear all
close all

id=figure;
hold on;
axis equal;

a=0.0;
b=1.0;
c=0.0;
d=1.0;
axis([a b c d]);

%grado spline
g=3;

%curva ancora non definita
init=0;

IM=1;
while (IM>0)
    IM=input('Vuoi 0=fine, 1=model, 2=deg-elev, 3=knot-ref, 4=load, 5=save ? ');
  
    if (IM==4)
        geom=read_nurbs('an_curve.txt');
        crv2=geom.nurbs;
        cla
        nrbctrlplot ( crv2 );
        qx=crv2.coefs(1,:);
        qy=crv2.coefs(2,:);
        mpk=length(qx);
        if (qx(1)==qx(mpk) & qy(1)==qy(mpk))
            chiusa=1;
        else
            chiusa=0;
        end
        nt=length(crv2.knots);
        g=crv2.order-1;
        if (crv2.knots(1)==crv2.knots(g+1) & crv2.knots(nt)==crv2.knots(nt-g))
            nagg=1;
        else
            nagg=0;
        end
        init=1;
    end
    
    if (IM==5 & init==1)
        nrbexport(crv2,'an_curve.txt');
    end
    
    if (IM==2 & nagg==1 & init==1)
      crv2 = nrbdegelev (crv2, 1);
      qx=crv2.coefs(1,:);
      qy=crv2.coefs(2,:);
      figure(id);
      nrbctrlplot ( crv2 );
      hold on;
      nrbkntplot( crv2 );
      axis equal
      axis([a,b,c,d]);
    end
    
    if (IM==3 & nagg==1 & init==1)
      crv2=nrbrefine(crv2,1,g-1);
      qx=crv2.coefs(1,:);
      qy=crv2.coefs(2,:);
      figure(id);
      nrbctrlplot ( crv2 );
      hold on;
      nrbkntplot( crv2 );
      axis equal
      axis([a,b,c,d]);
    end
%     if (nagg==0)
%       disp('Ki e De non disponibili per nodi aggiuntivi distinti')
%     end

    while (IM==1)
     if(init==0)
        [crv2,nagg,chiusa]=from_scratch(g);
        nrbctrlplot ( crv2 );
        g=crv2.order-1;
        qx=crv2.coefs(1,:);
        qy=crv2.coefs(2,:);
        init=1;
     else
      [u,v] = ginput(1);
      if (a<= u && u<=b & c<= v & v <=d)
       gc_ind=an_dist(u,v,qx,qy);
       [u,v] = ginput(1);
       plot (u,v,'+');
       for i=1:length(gc_ind)
         j=gc_ind(i);
         qx(j) = u;
         qy(j) = v;
       end
       crv2.coefs(1,:)=qx;
       crv2.coefs(2,:)=qy;
       figure(id);
       hold off;
       nrbctrlplot ( crv2 );
       axis equal
       axis([a,b,c,d]);
       hold on;
      else
         IM=4;
      end
     end
    end

end

function [an_ind,mind]=an_dist(x0,y0,x,y)
n=length(x);
mind=sqrt((x0-x(1)).^2+(y0-y(1)).^2);
m=1;
an_ind(m)=1;
for j=2:n
    tempd=sqrt((x0-x(j)).^2+(y0-y(j)).^2);
    if (tempd <= mind)
        if (tempd == mind)
            m=m+1;
        else
            m=1;
        end
        an_ind(m)=j;
        mind=tempd;
    end
end
end

function [crv2,nagg,chiusa]=from_scratch(g)
%Input interattivo dei punti di controllo della curva
%g --> grado spline
chiusa=0;
ncp=0;
button=1;
qx=[];
qy=[];
while button==1
   [xx,yy,button]=ginput(1);
   if (button==1) | (button==2)
    if button==2
        chiusa=1;
        xx=qx(1);
        yy=qy(1);
    end
    if (ncp>1)
        [gc_ind,dist]=an_dist(xx,yy,qx,qy);
        ncp=ncp+1;
        if(dist < 1.0e-2)
           qx(ncp)=qx(gc_ind(1));
           qy(ncp)=qy(gc_ind(1));
        else
           qx(ncp)=xx;
           qy(ncp)=yy;
        end
    else
        ncp=ncp+1;
        qx(ncp)=xx;
        qy(ncp)=yy;
    end
    plot(qx(ncp),qy(ncp),'ro');
   end
end 

%definizione struttura curva spline/nurbs
ctrl=zeros(4,ncp);
ctrl(1,:)=qx;
ctrl(2,:)=qy;
ctrl(3,:)=zeros(1,ncp);
ctrl(4,:)=ones(1,ncp);
% %grado spline
% g=3;
nt=ncp+g+1;
nagg=0;
for i=1:g
    if(qx(i)~=qx(ncp-g+i) | qy(i)~=qy(ncp-g+i))
        nagg=1;
    end
end
%nodi equispaziati e nodi aggiuntivi distinti
if (nagg==0)
    t=linspace(0,1,nt);
end
%nodi equispaziati e nodi aggiuntivi coincidenti
if (nagg==1)
    t(1:g)=zeros(1,g);
    t(g+1:nt-g)=linspace(0,1,nt-2*g);
    t(ncp+2:nt)=ones(1,g);
end

%definizione curva spline
crv2 = nrbmak (ctrl, t);
end
