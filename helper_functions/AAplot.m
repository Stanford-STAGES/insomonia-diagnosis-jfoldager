function projection = AAplot(XC,S,K,labels,labelNames)
figure
hb = bar(XC');
legend(labelNames{1},'Location','northeastoutside');
xlabel('Archetype')
xticklabels({'GS','Light sleeper','Long arousals', ...
    'Short sleeper','High ARI, AHI','Light sleep, PLMI ','Low TST, SE. High SOL'})
% 'GS','Low TST, SE','High ARI, AHI', ...
%     'Long arousals','LS, low SE, TST','Long SOL','High PLMI'
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)

xtickangle(90)
ylabel('Coordinate')
box off, grid on

deltaPhi    = 360/K;
orig        = [0,0];
idx         = 1;
figure
    hold on
    plot(orig(1), orig(2),'*'), axis equal
    clear cors
for k = 0:deltaPhi:360-deltaPhi
    cors(idx,:) = [sind(k),cosd(k)];
    txt1 = sprintf(' A$_%i $',idx);
    plot(orig(1), orig(2),sind(k),cosd(k),'*')
    text(sind(k),cosd(k),txt1)
    hold on
    idx = idx + 1;
end
hold on
plot([cors(:,1);cors(1,1)],[cors(:,2);cors(1,2)],'b')
hold on
projection = (S')*cors;
if ~isempty(labels) && isinteger(labels)
    gscatter(projection(:,1),projection(:,2),labels);
elseif ~isempty(labels) && ~isinteger(labels)
    scatter(projection(:,1),projection(:,2));
else
    scatter(projection(:,1),projection(:,2));
end

ax = gca;
ax.Visible = 'off';
figure; clf;

% Define proxies for data
x = projection(:,1);
y = projection(:,2);
z_orig = labels;
dyn_range_correction = min(z_orig)+1;
z = z_orig - dyn_range_correction;

% Define some arbitrary placement for archetypes and their names
XC = cors';% [0, 1; .875, -.5; -.875, -.5]'; % archetypes
atypenames = labelNames{2};
titlestr = 'Fitted plane';

noc = size(XC,2); % number of components
n_samples = 1e6; % number of points to generate surface from

% Generate samples in archetypal space
S = log(rand(noc, n_samples)); 
S = S ./ repmat(sum(S,1), [noc, 1]);
S = [S, eye(noc)];
XCS = XC*S;

% Fit surface to known data points in archetypal space for a given feature
sf = fit([x, y], z, 'poly23');

% Evaulate surface at known datapoints and determine triangulation for
% plotting
z_surf = feval(sf, [x, y]);
tri_data = delaunay(x, y);

% Or: use random samples in the archetypal space to be the basis for
% surface
z_surf_sampling = feval(sf, XCS');
tri_sampling = delaunay(XCS(1,:), XCS(2,:));

subplot(1,1,1)
    set(gca,'View', [55, 16])
    hold on
    % Plot data in archetypal "plane"
    scatter3(x, y, 0*z, 1, 'k');
    
    % Base surface on actual data points
    %trisurf(tri_data, x, y, z_surf, 'EdgeColor', 'none');
    
    % Or: base surface on sampling in archetypal space:
    
    trisurf(tri_sampling, XCS(1,:), XCS(2,:), z_surf_sampling', 'EdgeColor', 'none');
    
    % Add labels for archetypes
    for i = 1:noc
        text(XC(1, i), XC(2, i), 0, atypenames{i}, ...
                    'FontSize', 18, ...
        'FontWeight','bold')
    end
    
    % Add "axis" in archetypal space
    
    zmax = max(z_surf_sampling(:));
    zmin = min(z_surf_sampling(:));
    for i = 1:noc
        line([XC(1,i) XC(1,i)], [XC(2,i) XC(2,i)] , [0 zmax], 'Color', 'k')
        text([XC(1,i)], [XC(2,i)] , [zmax], num2str(round(zmax+dyn_range_correction,2)), ...
                    'FontWeight','bold')
        text([XC(1,i)], [XC(2,i)] , [zmin], num2str(round(zmin+dyn_range_correction,2)), ...
                    'FontWeight','bold')
        %scatter3([XC(1,i)], [XC(2,i)] , [zmin], 'ko')
        %scatter3([XC(1,i)], [XC(2,i)] , [zmax], 'ko')
        j = mod(i,noc)+1;
        line([XC(1,i) XC(1,j)], [XC(2,i) XC(2,j)] , [zmax zmax], 'Color', 'k')
        line([XC(1,i) XC(1,j)], [XC(2,i) XC(2,j)] , [0 0], 'Color', 'k')
    end
    
    title(titlestr)
    
    % Setup axis and limits
    axis equal
    axis off
ax = gca;
ax.Visible = 'off';

%%

end