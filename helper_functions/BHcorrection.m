function [corrected] = BHcorrection(pVals, Q)
    % EACH PHENOTYPE
    m   = size(pVals,1);
    i   = 1:m;
    pBH = (i./m)'.*Q;
    corrected = pVals.*repmat((m./i)',1,size(pVals,2));
    for idx = 1:size(pVals,2)
        pCompare(:,idx) = pVals(:,idx) < pBH';
    end
end