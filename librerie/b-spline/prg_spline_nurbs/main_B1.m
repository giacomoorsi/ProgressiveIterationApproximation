%script main_B1.m
clear all
close all
figure(1)
hold on;
ctrl=[0 2 4 6; 
      0 1 0 2; 
      0 0 0 0; 
      1 1 1 1];
crv2 = nrbmak (ctrl, [0 0 0 0.5 1 1 1]);
nrbctrlplot ( crv2 );
nrbkntplot( crv2 );
crv2
crv2.coefs
crv2.knots