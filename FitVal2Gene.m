function gene = FitVal2Gene( xInput )
gene = zeros(12, 1);
x = floor( xInput*1000*1000);
ret = mod(x, 10);
for i = 1:12
    gene(i) = ret;
    x = floor( (x-ret)/10 );
    ret = mod(x, 10);
end

end

