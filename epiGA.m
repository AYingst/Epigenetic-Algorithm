%   EpiGenetic Algorithm
%   Andrew Yingst
%	Purpose: expand the genetic algorithm to detect and account for
%	environmental influences  
close all; clf; clc; clear;

maxGens = 200; % maximum numbers of generations
bestNum = 50; % number of guaranteed most fit selected, must be even
popSize = 250; % num of individuals = size of popultion
checkGen = 6; % when to check for clusters
corrLimit = 0.7; % how high correlation must be
consistLim = 1.0; % inconsistency measurement threshold for clustering declaration

doClust = true; % initialize variable to look for clusters in the population indicating multiple optima
% Visualize workspace
chrom = 0:1/popSize:1-1/popSize; % chromosome distribution
chrom = repelem(10*chrom, size(chrom, 2), 1);% make square matrix for 3D plotting
envVar = 0:1/popSize:1-1/popSize; % environmental variable
envVar = (envVar/5 - .1)'; % resize and realign
envVar = repelem(envVar, 1, size(envVar, 1)); % make square matrix for 3D plotting
envOffset = chrom(1, round(size(chrom, 2)/2));% x-offset for symmetry; to pass to fitness func
envEff = envVar.*chrom - envVar*envOffset; % environmental effect
RandEnvEff = .2*rand(popSize, 1)-.1; % random environmental variable
fitMap = envEff + sin(chrom); % basic sine function + environmental effect
figure(1);
movegui('west');
s = surf(chrom, envVar, fitMap);
s.EdgeColor = 'none';
title('Chromosomal Fitness Effected by Environment', 'FontSize',16);
title('Environmental Impact');
xlabel('Chromosome Value');
ylabel('Environmental Variable');
zlabel("Fitness");
hold on;

%% Randomly generate an initial population of individuals
individual = zeros(popSize, 4);
for i = 1: 1: popSize
    individual(i,1) = i;% individual identifier
    individual(i,2) = 10*rand; % random chromosome 0:10
    individual(i,3) = .2*rand-.1; % random environmental variable -.1:.1
    % next 2 lines are initially unused epigenetic component of chromosome
    individual(i,4) = 0; %new cluster
    individual(i,5) = Inf; %activation boundary
end
%%
for numGens = 1:maxGens % number of generations for population reproduction
    %% calculate fitness of each individual
    fitness = zeros(popSize, 2);
    for i = 1:popSize
        fitness(i, 1) = i;
        % compare environment variable to boundary and modify chrom by cluster...
        % offset if necessary
        if (individual(i, 3) > individual(i, 5)) && (numGens >= checkGen + 1) % activate epigene if above threshold
            modChrom = individual(i, 4);
        else
            modChrom = individual(i, 2);
        end
        % calculate fitness for environmentally selected solution
        fitness(i, 2) = sin(modChrom) + individual(i,3)*modChrom...
            - individual(i,3)*envOffset;
    end
    %%
    figure(1); % redraw surface to erase past generations
    s = surf(chrom, envVar, fitMap);
    s.EdgeColor = 'none';
    view(-0.334,90);
    hold on;
    for i = 1:size(individual, 1)% plot each individual
        if (individual(i, 3) > individual(i, 5)) && (numGens >= checkGen +1) % activate epigene if above threshold
            modChrom(i,1) = individual(i, 4);
        else
            modChrom(i,1) = individual(i, 2);
        end
    end
    plot3(modChrom(:,1), individual(:,3), fitness(:,2), 'kx');
    title("Map of Genetic Solutions, Generation #" + numGens, 'FontSize',16);
    xlabel('Chromosome Value');
    ylabel('Environmental Variable');
    zlabel('Fitness');
    hold off;
    %%
    % detect clusters by tracking consistency of largest 2 clusters x6
    % generations. Ignore initial random population and stop after 5
    % samples
    if  (numGens < (checkGen + 1)) && (numGens >= 2)
        pause(.5); % for visual plot change analysis
        % calculate distance between all chromosomes
        chPairDist = pdist(individual(:,2));
        % detect clusters by linking chromosome groups into ever larger
        % groups
        chLinkList = linkage(chPairDist,'centroid');
        % measure the consistency of each cluster separation compared to
        % smaller clusters
        consistList = inconsistent(chLinkList);
        % the last entry is the consistency of the remaining 2 clusters.
        % low number indicates randomness forced into clusters
        if (consistList(end,4) < consistLim) 
            % do epigenetic modified GA only if real clustering exists for 5 generations
            doClust = false;
        end
        figure(2);
        movegui('east');
        dendrogram(chLinkList);
        title("Dendrogram of Generation #" + numGens, 'FontSize',16);
        xlabel('Cluster Identifier');
        ylabel('Cluster Separation');
        disp("Inconsistency " + consistList(end,4));
    end
    %%
    if (numGens == checkGen) && (doClust == true)
        % derive separation of clusters from consistency calculation
        chClustSep = chLinkList(end,3);
        chromCluster = cluster(chLinkList,'MaxClust',2);
        [bestChrome,I] =sort(fitness(:,2), 'descend');% sort by most fit
        % confirm clusters are related to environmental variables
        envCorr = corrcoef(chromCluster(I(1:bestNum),1), individual(I(1:bestNum),3));%correlation between most fit and envVar
        disp("envCorr = " + envCorr(1, 2));
        falseVar = .2*rand(bestNum, 1)-.1; % random environmental variable unrelated to fitness
        falseEnvCorr = corrcoef(bestChrome(1:bestNum),falseVar);%correlation between most fit and false envVar
        disp("falseEnvCorr = " + falseEnvCorr(1, 2));
        % convert all individuals from both clusters into uniform chromosome with
        % epigenetic programming
        if abs(envCorr(1, 2)) > corrLimit
            for i = 1:size(chromCluster, 1) % convert old chrom to new
                if individual(i,2) > chClustSep
                    individual(i, 4) = individual(i, 2);
                    individual(i, 2) = individual(i, 4) - chClustSep;
                else
                    individual(i, 4) = individual(i, 2) + chClustSep;
                end
            end
        end
        % determine environment boundary between 2 clusters and add mutation
        individual(:,5) = (max(individual(I(1:bestNum),3))+min(individual(I(1:bestNum),3)))/2 + .01*rand(popSize,1)-.005;
    end
    %%
    % Obtain the parent group (higher fitness = higher probability of being selected as the parent)
    [bestChrome,I] = sort(fitness(:,2), 'descend');% sort by most fit
    parent = zeros(popSize/2, 2);
    for i= 1:popSize/2
        if i <= bestNum/2 % select bestNum most fit individuals, rest come from roulette wheel
            parent(i,1) = I(2*i-1);
            parent(i,2) = I(2*i);
        else
            parent(i,:) = SelectParents(fitness);% Each line = the number of the two parents
        end
    end
    % Cross parents at random cut
    tempP = zeros(popSize/2, 2);
    for i = 1:popSize/2
        genesOfTwoChildren = CrossFit(parent(i,:), individual, checkGen, numGens);%return column vector = genes of two children
        geneOfChild1 = genesOfTwoChildren(1:12, 1:4);% 12 gene
        geneOfChild2 = genesOfTwoChildren(13:24, 1:4);% 12 gene
        tempP(2*i-1, 1) = 2*i-1;
        tempP(2*i-1, 2) = Gene2FitVal(geneOfChild1(:,1));% where the individual is located
        tempP(2*i-1, 4) = Gene2FitVal(geneOfChild1(:,2));% separation and boundary stay with child of most significant part parent gene
        tempP(2*i-1, 5) = Gene2FitVal(geneOfChild1(:,3));
        if geneOfChild1(1, 4) == 1
            tempP(2*i-1, 5) = -tempP(2*i-1, 5);
        end
        tempP(2*i,  1) = 2*i;
        tempP(2*i,  2) = Gene2FitVal(geneOfChild2(:,1));% where the individual is located
        tempP(2*i, 4) = Gene2FitVal(geneOfChild2(:,2));% separation and boundary stay with child of most significant part parent gene
        tempP(2*i, 5) = Gene2FitVal(geneOfChild2(:,3));
        if geneOfChild2(1, 4) == 1
            tempP(2*i, 5) = -tempP(2*i, 5);
        end
    end

    %Load new generation
    individual=tempP;
    individual(:,3) = .2*rand(popSize, 1)-.1; % random environmental variable added back in
    % data used for graphing progress
    data(numGens,:) = [min(individual(:,2)), max(individual(:,2)), mean(individual(:,2)),...
                       min(individual(:,5)), max(individual(:,5)), mean(individual(:,5)),...
                       min(individual(:,4)), max(individual(:,4)), mean(individual(:,4))];
%pause; %uncomment for slow progression
end
figure(3);

hold on;
errorbar(7:maxGens, data(7:maxGens, 3), data(7:maxGens, 3) - data(7:maxGens, 1),...
        data(7:maxGens, 2) - data(7:maxGens, 3) );
title("Refinement of Chromosome over Generations", 'FontSize',20);
        xlabel('Generation Number', 'FontSize',16);
        ylabel('Solution 1', 'FontSize',16);
hold off;
figure(4);
hold on;
errorbar(7:maxGens, data(7:maxGens, 6), data(7:maxGens, 6) - data(7:maxGens, 4),...
        data(7:maxGens, 5) - data(7:maxGens, 6) );
title("Refinement of Chromosome over Generations", 'FontSize',20);
        xlabel('Generation Number', 'FontSize',16);
        ylabel('Boundary Identifier', 'FontSize',16);
hold off;
figure(5);
hold on;
errorbar(7:maxGens, data(7:maxGens, 9), data(7:maxGens, 9) - data(7:maxGens, 7),...
        data(7:maxGens, 8) - data(7:maxGens, 9) );
title("Refinement of Chromosome over Generations", 'FontSize',20);
        xlabel('Generation Number', 'FontSize',16);
        ylabel('Solution 2', 'FontSize',16);
hold off;
