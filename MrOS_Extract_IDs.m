startup
extDiskPSGPath  = '';
extDiskAnoPath  = '';
listing         = extractfield(dir(extDiskPSGPath), 'name')';
edflisting      = listing(contains(lower(listing),'.edf'));
listing         = extractfield(dir(extDiskAnoPath), 'name')';
anolisting      = listing(contains(lower(listing),'.xml'));

% Find differences
% A = split(edflisting,'.');
% A = A(:,1);
% B = split(anolisting,'-');
% B = B(:,1:3);
% B = join(B,'-');
% C = find(contains(B,A)== 0);
% D = B(C);


% Synthesize data
A = split(anolisting,'-');
B = A(:,1:3);
BB = A(:,3);
C = join(B,'-');
s='A':'Z';
S = mat2cell(s(randi([1,length(s)],length(C),1)), 1, ones(length(C), 1));
D = [C, BB, num2cell(rand(length(C),1)), ...
    num2cell(randi([1,100],length(C),1)),...
    S'];
T = cell2table(D);
T.Properties.VariableNames = {'ID_Original','ID_Unique','Phenotype_1','Phenotype_2','Phenotype_3'};
writetable(T,'MrOS_Dataset.csv');