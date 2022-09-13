% Some more tests of the inpt forcing file type

%% Test 1: melt rate and no thickness update
%test whether the melt follows the input forcing
run = "1007"; %sinusoidal forcing
%run = "1008"; %stochastic forcing

fname = strcat("./ATTR_" ,run, "/outputs/outfile.nc"); %am
pyc = load(strcat("./ATTR_" ,run, "/pyc_centre.mat"));

%input file quantities
time = ncread(fname, 'TIME');
h = ncread(fname, 'h');
grfrac = ncread(fname, 'grfrac');
melt = ncread(fname, 'm');
grounded_volume = zeros(1,length(time));
mean_melt = zeros(1,length(time));
for i = 1:length(time)
    hnow = squeeze(h(:,:,i));
    grfracnow = squeeze(grfrac(:,:,i));
    grounded_volume(i)  = sum(sum(hnow .* grfracnow .* dx * dy));
    idx = grfracnow == 0;
    meltnow = squeeze(melt(:,:,i));
    mean_melt(i) = mean(meltnow(idx));
end

figure(1); clf; plot(time, mean_melt, 'ro-')
yyaxis right
plot(pyc.t, pyc.pc, 'bo-')
%% Look at the contours of melt
if 0
figure(2); clf
mmax = max(max(max(melt)));
mmin = min(min(min(melt))); %max and min for saturation
for i = 1:length(time)
    m = squeeze(melt(:,:,i))';
    m(1,1) = mmax;
    m(1,2) = mmin;
    clf; contourf(m, 30, 'linestyle', 'none');
    title(['t = ' num2str(time(i))])
    colorbar
    drawnow
    %pause
end
end