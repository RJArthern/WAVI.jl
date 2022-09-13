%Compare forcing types: input forcing file (where ambient temperature and
%salinity are specified as inpt files) and analytic qadratic forcing, wehre
%inpt specified via analytic temp and salinity profiles
%% Simple melt rate at time t = 0
fname1 = "./ATTR_1004/outputs/outfile.nc";
fname2 = "./ATTR_1003/outputs/outfile.nc";
m_newform = ncread(fname1, 'm');
m_newform = squeeze(m_newform(:,:,1));
m_oldform = ncread(fname2, 'm');
m_oldform = squeeze(m_oldform(:,:,1));
mdiff = m_newform - m_oldform;


figure(1); clf; subplot(3,1,1); contourf(m_newform', 30, 'linestyle', 'none'); colorbar; title('forcing file')
subplot(3,1,2); contourf(m_oldform', 30, 'linestyle', 'none'); colorbar; title('analytic forcing')
subplot(3,1,3); contourf(mdiff', 30, 'linestyle', 'none'); colorbar; title('difference (forcing file - analytic)')

%% 100 year run
fname_file = "./ATTR_1005/outputs/outfile.nc";
fname_analytic = "./ATTR_1006/outputs/outfile.nc";



%input file quantities
time_file = ncread(fname_file, 'TIME');
h_file = ncread(fname_file, 'h');
grfrac_file = ncread(fname_file, 'grfrac');
melt_file = ncread(fname_file, 'm');
grounded_volume_file = zeros(1,length(time_file));
mean_melt_file = zeros(1,length(time));
for i = 1:length(time_file)
    h = squeeze(h_file(:,:,i));
    grfrac = squeeze(grfrac_file(:,:,i));
    grounded_volume_file(i)  = sum(sum(h .* grfrac .* dx * dy));
    idx = grfrac == 0;
    melt = squeeze(melt_file(:,:,i));
    mean_melt_file(i) = mean(melt(idx));
end

%input file quantities
time_analytic = ncread(fname_analytic, 'TIME');
h_analytic = ncread(fname_analytic, 'h');
grfrac_analytic = ncread(fname_analytic, 'grfrac');
grounded_volume_analytic = zeros(1,length(time_analytic));
melt_analytic = ncread(fname_analytic, 'm');
mean_melt_analytic = zeros(1,length(time));
for i = 1:length(time_analytic)
    h = squeeze(h_analytic(:,:,i));
    grfrac = squeeze(grfrac_analytic(:,:,i));
    grounded_volume_analytic(i)  = sum(sum(h .* grfrac .* dx * dy));
    idx = grfrac == 0;
    melt = squeeze(melt_analytic(:,:,i));
    mean_melt_analytic(i) = mean(melt(idx));
end

figure(2); clf;
subplot(2,2,1);
plot(time_analytic, mean_melt_analytic, 'r', time_file, mean_melt_file, 'bo')
legend('analytic', 'file');
title('Mean melt rate');

subplot(2,2,2);
plot(time_analytic, (mean_melt_analytic - mean_melt_file)./mean_melt_file * 100, 'ko-');
title('Mean melt rate % diff');

subplot(2,2,3);
plot(time_analytic, grounded_volume_analytic/1e12, 'r', time_file, grounded_volume_analytic/1e12, 'bo')
legend('analytic', 'file');

subplot(2,2,4);
plot(time_analytic, (grounded_volume_analytic - grounded_volume_file)./grounded_volume_file * 100, 'ko-');
title('grounded volume % diff');