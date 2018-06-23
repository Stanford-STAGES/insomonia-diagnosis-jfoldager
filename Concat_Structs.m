startup
extDiskPSGPath = '';
a = load(strcat(extDiskPSGPath,'features'));
b = load(strcat(extDiskPSGPath,'One2_1099_530'));
f = fieldnames(b.FEATURES);
for i = 1:length(f)
    a.FEATURES.(f{i}) = b.FEATURES.(f{i});
end
MrOs_Features = a.FEATURES;
MrOs_Features = orderfields(MrOs_Features);
save(strcat(''),'MrOs_Features')