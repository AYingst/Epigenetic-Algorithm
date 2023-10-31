% Convert fitness value to gene then cross parent genes at random cut position
function childGene = CrossFit(parent, individual, checkGen, numGens)
childGene = zeros(24, 4);
parent1Index = parent(1);
parent2Index = parent(2);
% select genome based on env
parent(1, 1) = individual(parent1Index, 2);
parent(2, 1) = individual(parent2Index, 2);
parent(1, 2) = individual(parent1Index, 4);
parent(2, 2) = individual(parent2Index, 4);
parent(1, 3) = individual(parent1Index, 5);
if individual(parent1Index, 5) < 0 % track and fix negative env variable
    parent(1, 3) = abs(parent(1, 3));
    childGene(4, 1) = 1;
end
parent(2, 3) = individual(parent2Index, 5);
if individual(parent2Index, 5) < 0
    parent(2, 3) = abs(parent(2, 3));
    childGene(4, 13) = 1;
end
% do crossover for 1 solution if GA or 3 if EGA
if numGens >= checkGen
    numOfCross = 3;
else
    numOfCross = 1;
end
for p = 1:numOfCross
    parent1Gene = FitVal2Gene(parent(1, p));
    parent2Gene = FitVal2Gene(parent(2, p));
    cut = floor( rand(1,1) * 12 ) + 1;
    for i=1:cut
        tmp = parent1Gene(i);
        parent1Gene(i) = parent2Gene(i);
        parent2Gene(i) = tmp;
        childGene(:, p) = [parent1Gene; parent2Gene];
    end

end
end
