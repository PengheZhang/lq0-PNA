function f = multi_logistic_fun(XW,Y)
[n,~] = size(Y);

if max( max( XW ) ) >= 10
    m_XW = max( XW,[],2 );
    XW = [ XW - m_XW, -m_XW ];
else
    XW = [ XW, zeros(n,1) ];
end
eXW = exp(XW);
allsum = sum( eXW,2 ); 

YXW = Y.*(XW);
f = ( sum(log(allsum)) - sum(sum(YXW)) )/n ;

end