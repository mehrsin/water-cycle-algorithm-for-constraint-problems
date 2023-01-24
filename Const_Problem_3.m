function c=Const_Problem_3(x)

nvars=10;

b=0;
for i=1:nvars
    bb=x(i)^2;
    b=bb+b;
end

c=b-1-eps;

end

