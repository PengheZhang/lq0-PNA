function [g,PW] = grad_Pw(Xt,Y,XW,n)
[~,m] = size(Y);
if max( max( XW ) ) >= 5
    m_XW = max( XW,[],2 );
    XW = [ XW - m_XW, -m_XW ];
else
    XW = [ XW, zeros(n,1) ];
end
eXW = exp(XW);
allsum = sum( eXW,2 ); 
PW = eXW( :, 1:m )./ allsum;  
g = Xt*( PW - Y )/n;
end