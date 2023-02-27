% Demonstration of lateral speckle tracking and channel correlation for
% estimation of transducer motion.
%
% Both methods here use a weighted least squares approach to multi-lag
% estimation over different frame pairs, using the cross-correlation 
% coefficient as a weight. Subsample estimation is used in both methods,
% here chosen as iterative reconstructive  
%
% Requires planewave_data.mat from the github release
%

%% Common parameters for estimation
clear all

depth = 20e-3; % Kernel depth (m)
k_axial = 10e-3; % Axial kernel length (m)
k_frame = 10; % Number of frames for multi-lag estimation

load('planewave_data.mat');

%% Channel correlation estimates
% Use raw RF channel data to estimate the lateral offset between frames 
% based on the observed shift in array subapertures

k_elem = 10; % Maximum number of element shifts to search in each direction

figure(1);clf
est_ch=lateral_channel_track_wls(rf,acq_params,...
    'depth',depth,'k_axial',k_axial,'k_elem',k_elem,'k_frame',k_frame,'make_plot',1);   
subplot(121)
hold all
plot(1:length(acq_params.locations),acq_params.locations,'k-')
axis tight
legend('Estimate','Truth')
hold off
title('Channel correlation')

print -dpng channel_correlation

%% Speckle tracking estimates

% Beamform images for use in speckle tracking
lambda=acq_params.c/acq_params.f0;
x=-10e-3:lambda/4:10e-3;
z=(depth-k_axial/2-1e-3):lambda/6:(depth+k_axial/2+1e-3);

[x1,z1]=meshgrid(x,z);
x1=x1(:);
z1=z1(:);

tx_dists=z1; % 0 degree PW
rx_dists=zeros(numel(z1),size(acq_params.rx_pos,1),'single');
for i=1:size(acq_params.rx_pos,1)
    rx_dists(:,i)=sqrt((acq_params.rx_pos(i,1)-x1).^2+...
                       (acq_params.rx_pos(i,3)-z1).^2);
end

% Focus by interpolation
t=(acq_params.t0+(0:size(rf,1)-1))/acq_params.fs;
data_focused=zeros(numel(x1),size(rf,3),'like',rf);
% NOTE: Change this to a parfor loop for faster computation
for i=1:size(rf,3) % Loop over transmissions
    for j=1:size(acq_params.rx_pos,1) % Loop over receive elements
        t_interp=(rx_dists(:,j)+tx_dists)/acq_params.c;
        data_focused(:,i)=data_focused(:,i)+interp1(t,rf(:,j,i),t_interp,'linear',0);
    end
end
data_focused=reshape(data_focused,length(z),length(x),size(rf,3));

figure(2);clf
k_lateral=10e-3;
k_search=2.5e-3;
env=abs(hilbert(data_focused));
est_sp=lateral_speckle_track_wls(env,z,x,...
    'depth',depth,'k_axial',k_axial,'k_lateral',k_lateral,'k_search',k_search,'k_frame',k_frame,'make_plot',1);   
subplot(121)
hold all
plot(1:length(acq_params.locations),acq_params.locations,'k-')
axis tight
legend('Estimate','Truth')
hold off
title('Speckle tracking')

print -dpng speckle_tracking