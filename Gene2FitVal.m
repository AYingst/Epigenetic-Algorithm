% Convert gene back to fitness value
function num = Gene2FitVal( gene )
num = 0;
for i = 1:12
	num = num + gene(i)*10^(i-7);
end
if num<0 || num >10
    num = 10*rand(1, 1);
end
end

