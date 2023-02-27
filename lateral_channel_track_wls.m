% est_avg = lateral_channel_track_wls(rf,acq_params,varargin)
%
% Use channel cross-correlation of kernels to perform lateral displacement estimation
% combine estimates made across different pairs of frames using a weighted
% least squares strategy
%
% Required:
% rf - ax x lat x frame
% acq_params - Structure containing fs, t0, c, rx_pos
% Name/value pair options:
% k_axial - Axial correlation kernel in m
% k_frame - Frames from reference to search
% k_elem - Maximum element shift to search
% depth - Depth to track at in m
% make_plot - Plot the estimates
%
function est_avg = lateral_channel_track_wls(rf,acq_params,varargin)

p=inputParser;
p.addParameter('k_axial',10e-3);
p.addParameter('k_frame',10);
p.addParameter('k_elem',10);
p.addParameter('depth',20e-3);
p.addParameter('make_plot',0);
p.parse(varargin{:});

rf=squeeze(single(rf));

if(p.Results.make_plot)
    clf
end

ksz=p.Results.k_axial/acq_params.c*2*acq_params.fs; %samples
depth=p.Results.depth/acq_params.c*2*acq_params.fs+acq_params.t0; %samples
n=round((-ksz/2:ksz/2)+depth);

% These two need to work together to define a maximum velocity
N=p.Results.k_frame; % How many frames to examine from a single reference
el=-p.Results.k_elem:p.Results.k_elem; % How many elements separations to search

% Define sampling matrix
H=zeros(1+size(rf,3)*N-N*(N+1)/2,size(rf,3));
inds=zeros(size(H,1),2);
H(1,:)=1/size(rf,3); % Zero mean condition
cnt=2;
for i=1:size(rf,3)
    for j=i+1:min(i+N,size(rf,3))
        H(cnt,i)=-1;
        H(cnt,j)=1;
        inds(cnt,:)=[i,j];
        cnt=cnt+1;
    end
end

% Make observation pairs
Y=zeros(size(H,1),1);
W=zeros(size(H,1),1);
Y(1)=0;
W(1)=1;
aperture=el*diff(acq_params.rx_pos([1 2],1));
dx=diff(aperture(1:2));
for i=2:size(inds,1)
    cc=zeros(1,length(el));
    for e=1:length(el)
        apods=false(2,size(rf,2));
        if(el(e)<0)
            apods(1,1:end+el(e))=1;
            apods(2,-el(e)+1:end)=1;
        else
            apods(2,1:end-el(e))=1;
            apods(1,el(e)+1:end)=1;
        end
        ref=rf(n,apods(1,:),inds(i,1));
        cur=rf(n,apods(2,:),inds(i,2));
        cc(e)=ncc(cur,ref);
    end
    
    [lag,rho]=subsample_peak(cc,el,'reconstructive_iterative');
    Y(i)=lag*dx;
    W(i)=rho; % can apply atanh(rho) here
end

% Weighted least squares
Wd=diag(W);
est_avg=(H'*Wd*H)\H'*Wd*Y;

if(p.Results.make_plot)
    subplot(121)
    plot(est_avg,'Color','r','LineWidth',2)
    xlabel('Frame')
    ylabel('Displacement (m)')
    subplot(122)
    histogram(W,'Normalization','probability')
    xlabel('Correlation')
    ylabel('Probability')
end
