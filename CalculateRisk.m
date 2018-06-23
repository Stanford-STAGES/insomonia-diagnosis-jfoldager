startup
load JansenDirections.mat
% Calculates risk based on effect directions from Jansen et al. (2018)
% Input: Insonia risk allele data, effect directions from Jansen
% Output: Sum over risk alleles to make up polygenic risk average
%% V1
load V1newdata.mat
for i = 1:size(X,2)
    if dircs(i) == 0
        X{:,i} = 2 - X{:,i};
    end
end
risks   = sum(X{:,:},2)/(2*size(X,2));
Y.Risks = risks;
figure, histogram(Y.Risks,'Normalization','probability');
xlabel('Polygenic Risk Coefficient'), ylabel('Probability')
grid on, box off
save('RisksV1','risks')
%% V2
load V2newdata.mat
load JansenDirections.mat
for i = 1:size(X,2)
    if dircs(i) == 0
        X{:,i} = 2 - X{:,i};
    end
end
risks   = sum(X{:,:},2)/(2*size(X,2));
Y.Risks = risks;
figure, histogram(Y.Risks);
save('RisksV2','risks')
