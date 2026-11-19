function main(inputFile,ALPHA,BETA0,RHO0,isINIT,INT,max_iters,iters_for_R,SEED)

% inputFile: the name of folder in 'Inputs/'. For example, 'DATA_02282025'.

% ALPHA: the value of \alpha in the main manuscript. For example, '2.0'.
% ALPHA should be of K/2 for a natural number K.

% BETA0: the initial value of \beta in the main manuscript.

% RHO0: the initial value of \rho in the main manuscript.

% isINIT: the training algorithm is designed to regularly store the interim
% results in the folder 'Initializations'. If it is 'Yes', then the
% algorithm will load the stored file and resume the training, and the 
% previously specified BETA0 and RHO0 will be ignored. If it is
% 'No', then the algorithm estimates the parameters from the very
% beginning, with BETA0 and RHO0 as the initialization.

% INT: the length of unitless interval. For example, '0.5'.

% max_iters: the maximum number of EM iterations. For example, '10000'.

% iters_for_R: the number of iterations between updates of the
% record-specific standardization parameter. If you want to fix the
% standardization parameter as its initialization in 
% 'Inputs/inputFile/average_sedrates.txt', set it to be 'inf'.

% SEED: random seed. If you want to "randomize" the running, set SEED to be
% 'NaN'.

if strcmp(isINIT,'Yes')
    INIT = [inputFile,'_K_',num2str(2.*ALPHA),'_INT_',num2str(INT)];
elseif strcmp(isINIT,'No')
    INIT = [];
end

if isnan(SEED) || strcmp(SEED,'NaN')
    SEED = 'shuffle';
end

AccMODEL(inputFile,ALPHA,BETA0,RHO0,INIT,INT,max_iters,iters_for_R,SEED);

end