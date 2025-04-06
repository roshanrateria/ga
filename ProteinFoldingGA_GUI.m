function ProteinFoldingGA_GUI
    %% Create the main GUI window
    f = figure('Position',[100,100,800,500],...
               'Name','GA for Protein Structure Prediction',...
               'NumberTitle','off',...
               'Resize','off');

    %% Panel for GA Parameter Inputs
    paramPanel = uipanel('Parent', f, 'Title', 'GA Parameters',...
                         'Position', [0.02 0.55 0.3 0.42]);
                     
    % Population size
    uicontrol('Parent', paramPanel, 'Style', 'text',...
              'Position', [10, 130, 120, 20],...
              'String', 'Population Size:','HorizontalAlignment','left');
    popSizeEdit = uicontrol('Parent', paramPanel, 'Style', 'edit',...
              'Position', [140, 130, 60, 25],...
              'String', '20');

    % Number of genes (angles)
    uicontrol('Parent', paramPanel, 'Style', 'text',...
              'Position', [10, 95, 120, 20],...
              'String', 'Number of Genes:','HorizontalAlignment','left');
    nGenesEdit = uicontrol('Parent', paramPanel, 'Style', 'edit',...
              'Position', [140, 95, 60, 25],...
              'String', '10');

    % Number of generations
    uicontrol('Parent', paramPanel, 'Style', 'text',...
              'Position', [10, 60, 120, 20],...
              'String', 'Generations:','HorizontalAlignment','left');
    maxGenEdit = uicontrol('Parent', paramPanel, 'Style', 'edit',...
              'Position', [140, 60, 60, 25],...
              'String', '50');

    % Mutation rate
    uicontrol('Parent', paramPanel, 'Style', 'text',...
              'Position', [10, 25, 120, 20],...
              'String', 'Mutation Rate:','HorizontalAlignment','left');
    mutRateEdit = uicontrol('Parent', paramPanel, 'Style', 'edit',...
              'Position', [140, 25, 60, 25],...
              'String', '0.1');

    % Crossover rate
    uicontrol('Parent', paramPanel, 'Style', 'text',...
              'Position', [10, -10, 120, 20],...
              'String', 'Crossover Rate:','HorizontalAlignment','left');
    crossRateEdit = uicontrol('Parent', paramPanel, 'Style', 'edit',...
              'Position', [140, -10, 60, 25],...
              'String', '0.8');

    %% Button to Start GA
    startBtn = uicontrol('Style','pushbutton',...
                         'String','Start GA',...
                         'FontSize',12,...
                         'Position',[50, 50, 200, 40],...
                         'Callback',@startGA);

    %% Axes for GA Fitness Progress
    fitnessAx = axes('Parent', f, 'Units','normalized',...
                     'Position',[0.35 0.55 0.6 0.42]);
    title(fitnessAx, 'GA Fitness Evolution');
    xlabel(fitnessAx, 'Generation');
    ylabel(fitnessAx, 'Best Fitness');

    %% Axes for Protein Structure Visualization
    structureAx = axes('Parent', f, 'Units','normalized',...
                     'Position',[0.1 0.05 0.8 0.35]);
    title(structureAx, 'Protein Structure Evolution');
    xlabel(structureAx, 'X');
    ylabel(structureAx, 'Y');
    grid(structureAx, 'on');
    
    %% Text field to show current best fitness
    bestFitnessTxt = uicontrol('Style','text',...
                               'FontSize',12,...
                               'Position',[300, 20, 200, 30],...
                               'String','Best fitness: ');

    %% Callback function for GA start button
    function startGA(~, ~)
         % Retrieve GA parameters from UI
         popSize = str2double(get(popSizeEdit, 'String'));
         nGenes = str2double(get(nGenesEdit, 'String'));
         maxGen = str2double(get(maxGenEdit, 'String'));
         mutationRate = str2double(get(mutRateEdit, 'String'));
         crossoverRate = str2double(get(crossRateEdit, 'String'));
         
         % Update structure axes limits based on number of genes (for proper scaling)
         axis(structureAx, [-nGenes-1, nGenes+1, -nGenes-1, nGenes+1]);
         
         % Initialize population: each individual is a vector of angles between -pi and pi
         pop = (rand(popSize, nGenes) * 2*pi - pi);
         bestFitness = inf;
         bestIndividual = [];
         fitnessHistory = zeros(maxGen,1);
         
         % Clear previous plots
         cla(fitnessAx);
         cla(structureAx);
         
         % GA main loop
         for gen = 1:maxGen
              % Evaluate fitness for each individual
              fitness = zeros(popSize,1);
              for i = 1:popSize
                  fitness(i) = fitnessFunc(pop(i,:));
              end

              % Keep track of the best solution in current generation
              [currentBest, idx] = min(fitness);
              if currentBest < bestFitness
                  bestFitness = currentBest;
                  bestIndividual = pop(idx,:);
              end
              fitnessHistory(gen) = bestFitness;
              
              % Update fitness plot every generation
              plot(fitnessAx, 1:gen, fitnessHistory(1:gen), '-o');
              xlabel(fitnessAx, 'Generation');
              ylabel(fitnessAx, 'Best Fitness');
              title(fitnessAx, 'GA Fitness Evolution');
              drawnow;

              % Visualize the current best protein structure
              coords = computeCoordinates(bestIndividual);
              plot(structureAx, coords(:,1), coords(:,2), '-o','LineWidth',2);
              xlabel(structureAx, 'X');
              ylabel(structureAx, 'Y');
              title(structureAx, sprintf('Protein Structure at Generation %d', gen));
              grid(structureAx, 'on');
              drawnow;
              
              % Update best fitness text
              bestFitnessTxt.String = sprintf('Best fitness: %.4f', bestFitness);
              
              % --- Selection using tournament selection ---
              newPop = zeros(size(pop));
              for i = 1:popSize
                   idxs = randi(popSize, [1,2]);
                   if fitness(idxs(1)) < fitness(idxs(2))
                        winner = pop(idxs(1),:);
                   else
                        winner = pop(idxs(2),:);
                   end
                   newPop(i,:) = winner;
              end

              % --- Crossover ---
              for i = 1:2:popSize-1
                   if rand < crossoverRate
                        cp = randi(nGenes-1); % Random crossover point
                        temp = newPop(i, cp+1:end);
                        newPop(i, cp+1:end) = newPop(i+1, cp+1:end);
                        newPop(i+1, cp+1:end) = temp;
                   end
              end

              % --- Mutation ---
              for i = 1:popSize
                   for j = 1:nGenes
                        if rand < mutationRate
                             newPop(i,j) = newPop(i,j) + randn*0.1;
                             % Keep gene within [-pi, pi]
                             newPop(i,j) = max(min(newPop(i,j), pi), -pi);
                        end
                   end
              end

              % Update population for next generation
              pop = newPop;
              
              % Delay to see the evolution (adjust pause time as desired)
              pause(0.2);
         end

         % Final plots update
         plot(fitnessAx, 1:maxGen, fitnessHistory, '-o');
         xlabel(fitnessAx, 'Generation');
         ylabel(fitnessAx, 'Best Fitness');
         title(fitnessAx, 'GA Fitness Evolution');
         
         coords = computeCoordinates(bestIndividual);
         plot(structureAx, coords(:,1), coords(:,2), '-o','LineWidth',2);
         xlabel(structureAx, 'X');
         ylabel(structureAx, 'Y');
         title(structureAx, 'Final Protein Structure');
         grid(structureAx, 'on');
         bestFitnessTxt.String = sprintf('Final best fitness: %.4f', bestFitness);
    end

    %% Fitness Function (minimize sum of squared angles)
    function fVal = fitnessFunc(individual)
         fVal = sum(individual.^2);
    end

    %% Compute 2D Coordinates for the Protein Chain Based on Angles
    % Starting at (0,0) and initial direction along the positive x-axis,
    % each angle rotates the direction by that amount.
    function coords = computeCoordinates(angles)
         n = length(angles);
         coords = zeros(n+1, 2);
         currentAngle = 0;
         % Starting point at origin (0,0)
         coords(1,:) = [0, 0];
         stepLength = 1; % constant step length for each segment
         for k = 1:n
              % Update direction by adding the current angle
              currentAngle = currentAngle + angles(k);
              % Compute new point
              newX = coords(k,1) + stepLength*cos(currentAngle);
              newY = coords(k,2) + stepLength*sin(currentAngle);
              coords(k+1,:) = [newX, newY];
         end
    end
end
