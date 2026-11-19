inputFile = 'EXPERIMENT';

% random seed:
RAND_SEED = 42; % if you want to randomize it, set RAND_SEED = NaN;

outputFile = Aligner(inputFile,RAND_SEED);

getFigures(outputFile);


