function [x,y] = epitrochoid(t)
% definita tra 0 e 6*pi
a=2;
b=3;
c=1;
x = (a-b)*cos(t)-c*cos((a/b+1)*t);
y = (a-b)*sin(t)-c*sin((a/b+1)*t);
end
