function f= Problem_3(x)

nvars=10;

a=1;
for i=1:nvars
    aa=x(i);
    a=aa*a;
end

f=-((sqrt(nvars))^nvars)*a;

end

