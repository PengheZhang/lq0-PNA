function Ay = Matvecn(Ainput,y)

yr = reshape(y,Ainput.nT,Ainput.m);
Xyr = Ainput.XT2*yr;
XTB = [];

for i = 1:Ainput.m
    B = Xyr(:,i).*Ainput.U{i,i};
    for j = i+1:Ainput.m
        B = B+Xyr(:,j).*Ainput.U{i,j};
    end
    for k = 1:i-1
        B = B+Xyr(:,k).*Ainput.U{k,i};
    end
   XTB = [XTB;Ainput.XT1*B];
end

Ay=XTB+Ainput.gam*y;

end