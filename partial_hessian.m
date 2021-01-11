function H = partial_hessian( PW, m_bar, n )

H=cell( m_bar, m_bar );

for i=1:m_bar
    H{i,i} = PW(:,i).*(1-PW(:,i))/n;
    for j = i+1:m_bar
        H{i,j} = - PW(:,i).*PW(:,j)/n;
    end
end
end