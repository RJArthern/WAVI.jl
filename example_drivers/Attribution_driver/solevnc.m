%Plot the evolution of a solution when files are listed as an ncfile

% preliminaries
folder = "./ATTR_1003/outputs/";
fname = strcat(folder,'outfile.nc');

%% grid info
dx = 4000;
dy = 4000;

time = ncread(fname, 'TIME');
h_stack = ncread(fname, 'h');
x = ncread(fname, 'x');
y = ncread(fname, 'y');
b = ncread(fname, 'b');
grfrac_stack = ncread(fname, 'grfrac');
melt = ncread(fname, 'm');

nx = length(x);
ny = length(y);
bathy = squeeze(b(:,:,1));

grounded_volume = zeros(1,length(time));
mean_melt = zeros(1,length(time));
for i = 1:length(time)
    h = squeeze(h_stack(:,:,i));
    grfrac = squeeze(grfrac_stack(:,:,i));
    grounded_volume(i)  = sum(sum(h .* grfrac .* dx * dy));
    idx = grfrac == 0;
    mm = sqeeze(melt(:,:,i));
    mean_melt(i) = mean(mm(idx));
end


%% Plots
% Plot 1: grounded fraction
figure(1); clf; 
plot(time, grounded_volume/1e9, 'ro-')
xlabel('time (years)');
ylabel('grounded volume (km^3)')
title('evolution of grounded volume');
%
%% Plot 2: evolution of the gl
figure(2); clf; hold on; box on
%  add base layer of bathymetry
contourf(x/1e3,y/1e3, bathy', 30, 'linestyle', 'none');
ax = gca;
colormap(ax, autumn(100)); 
xlabel('x (km)');
ylabel('y (km)');
title('bathymetry and grounding line position')


% add gronding lines
axnew = axes;
axnew.Position = ax.Position;
hold on

colmap = parula(length(time));
for i = 1:length(time)
    contour(x/1e3,y/1e3,squeeze(grfrac_stack(:,:,i))', [0.5, 0.5], 'linecolor', colmap(i,:));
   % drawnow
   % pause
end

c = colorbar(axnew);
c.TickLabels = {'0', num2str(time(end))};
c.Ticks = [min(c.Ticks), max(c.Ticks)];
c.Label.String = 'time (years)';
axnew.Position = ax.Position;
axnew.Visible = 'off';

%% Plot 3: slice along the centre
figure(3); clf; hold on; box on
idx = floor(ny/2);
b = bathy(:,idx);
plot(x, b, 'm', 'linewidth', 2);

for i = 1:length(time)
    h = squeeze(h_stack(:,idx,i));
    grfrac = squeeze(grfrac_stack(:,idx,i));
    isfloat = (grfrac == 0);
    
    base = zeros(size(x));
    base(~isfloat) = b(~isfloat);
    base(isfloat) = -918/1028 * h(isfloat);
    s = h + base;
    plot(x,s, 'color', colmap(i,:));
    plot(x,base, '--', 'color', colmap(i,:));
end
c = colorbar;
c.TickLabels = {'0', num2str(time(end)/2), num2str(time(end))};
c.Ticks = [min(c.Ticks), (min(c.Ticks) + max(c.Ticks))/2,  max(c.Ticks)];
c.Label.String = 'time (years)';
xlim([min(x), max(x)]);
