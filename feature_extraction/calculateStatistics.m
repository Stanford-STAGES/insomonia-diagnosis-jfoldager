function out = calculateStatistics(in,hypnogram,qs)
    ss                  = [0,1,2,3,5];
    ssn                 = {'WA';'N1';'N2';'N3';'RE'};
    qs                  = [0.25, 0.50, 0.75];
    qn                  = strcat(repmat(in.Properties.VariableNames,size(qs,2),1),...
        repmat({'q1';'q2';'q3'},1,size(in,2)));
    qn                  = qn(:)';
    out                 = table();
    for s = 1:length(ss)
        qv = quantile(in{hypnogram == ss(s),:},qs);
        vals = qv(:)';
        names = strcat(qn,repmat(ssn(s),1,size(qn,2)));
        out     = [out, array2table(vals,'VariableNames',names)];
    end
end