function crv=nrbrefine(crv,nctrlptsins,regularity)
[~,~,new_knots]=kntrefine(crv.knots,nctrlptsins,crv.order-1,regularity);
% disp('ik')
% disp(new_knots)
crv = nrbkntins (crv , new_knots);
end