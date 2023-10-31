% Roulette wheel selection of most fit individuals
% Duplication is allowed
function parents = SelectParents( fitness )
%adjust fitness scores to all positive
parents = zeros(1, 2);
% account for negative numbers
minFit = min(fitness);
minFit = minFit(:,2);
if minFit < 0 %shift range up
    fitness(:,2) = fitness(:,2) - minFit;
end 
sumfit = 0;
for i = 1:size(fitness, 1)
    sumfit = sumfit + fitness(i,2);
end
randomValue = rand(1, 2);
possibility = randomValue .* sumfit;
for parentNum = 1:2
    tmpSum = 0;
    for i = 1:size(fitness, 1)
        tmpSum = tmpSum + fitness(i,2);
        if tmpSum >= possibility(parentNum)
            parents(parentNum) = fitness(i,1);
            break;
        end
    end 
end
end
