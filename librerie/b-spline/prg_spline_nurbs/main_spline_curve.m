%Questo programma permette di modellare una curva spline 2D
% si danno i punti interattivamente con button1; dato l'ultimo punto
% si prema button3 per fine inserimento e curva aperta, oppure button2
% per fine inserimento e curva chiusa
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

%Input interattivo dei punti di controllo della curva
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

%definizione struttura curva spline
ctrl=zeros(4,ncp);
ctrl(1,:)=qx;
ctrl(2,:)=qy;
ctrl(3,:)=zeros(1,ncp);
ctrl(4,:)=ones(1,ncp);

%grado spline
g=3;

%numero nodi partizione estesa
nt=ncp+g+1;
% t=linspace(0,1,nt);

%definizione partzione estesa a noid aggiuntivi coincidenti
t(1:g)=zeros(1,g);
t(g+1:nt-g)=linspace(0,1,nt-2*g);
t(ncp+2:nt)=ones(1,g);

%definizione curva spline
crv2 = nrbmak (ctrl, t);
nrbctrlplot ( crv2 );

IM=1;
while (IM>0)
    IM=input('Vuoi 0=uscire, 1=modellare ? ');

    while (IM==1)
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
         IM=2;
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
