%script main_B4.m
clear all
close all
figure(1)
hold on;

crv = nrbcirc (2, [0 0 0], 0, 2*pi); 
nrbctrlplot ( crv );
nrbplot( crv, 55 );
crv.coefs
crv.knots

%TO DO