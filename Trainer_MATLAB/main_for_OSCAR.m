function main_for_OSCAR(inputFile,ALPHA,BETA0,RHO0,isINIT,INT,max_iters,iters_for_R)

if strcmp(isINIT,'Yes')
    INIT = [inputFile,'_K_',num2str(2.*ALPHA),'_INT_',num2str(INT)];
elseif strcmp(isINIT,'No')
    INIT = [];
end

AccMODEL(inputFile,ALPHA,BETA0,RHO0,INIT,INT,max_iters,iters_for_R);

end