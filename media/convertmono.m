function Yout = convertmono (Y)
if size(Y,2) == 2
    Yout = (Y(:,1) + Y(:,2))./2;
else
    Yout = Y;
end
