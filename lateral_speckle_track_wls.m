% est_avg = lateral_speckle_track_wls(env,z,x,varargin)
%
% Use lateral cross-correlation of kernels to perform speckle tracking and
% combine estimates made across different pairs of frames using a weighted
% least squares strategy
%
% Required:
% env - ax x lat x frame
% z - Axial distance vector (m)
% x - Lateral distance vector (m)
% Name/value pair options:
% k_axial - Axial correlation kernel in m
% k_frame - Frames from reference to search
% k_lateral - Lateral kernel size in m
% k_search - Lateral distance to search in m
% depth - Depth to track at in m
% make_plot - Plot the estimates
%
function est_avg = lateral_speckle_track_wls(env,z,x,varargin)

p=inputParser;
p.addParameter('k_axial',10e-3);
p.addParameter('k_lateral',10e-3);
p.addParameter('k_frame',10);
p.addParameter('k_search',5e-3);
p.addParameter('depth',20e-3);
p.addParameter('make_plot',0);
p.parse(varargin{:});

env=squeeze(single(env));

if(p.Results.make_plot)
    clf
end

assert(p.Results.k_lateral+p.Results.k_search*2<=diff(x([1 end])),'Lateral extent too small for kernel and search');
assert(p.Results.k_axial<=diff(z([1 end])),'Axial extent too small for kernel');
N=p.Results.k_frame; % How many frames to examine from a single reference
xinds=find(abs(x-mean(x))<p.Results.k_lateral/2);
zinds=find(abs(z-p.Results.depth)<p.Results.k_axial/2);
n=ceil(p.Results.k_search/diff(x(1:2)));
srch=(-n:n);

% Define sampling matrix
H=zeros(1+size(env,3)*N-N*(N+1)/2,size(env,3));
inds=zeros(size(H,1),2);
H(1,:)=1/size(env,3); % Zero mean condition
cnt=2;
for i=1:size(env,3)
    for j=i+1:min(i+N,size(env,3))
        H(cnt,i)=1;
        H(cnt,j)=-1;
        inds(cnt,:)=[i,j];
        cnt=cnt+1;
    end
end

% Make observation pairs
Y=zeros(size(H,1),1);
W=zeros(size(H,1),1);
Y(1)=0;
W(1)=1;
dx=diff(x(1:2));
for i=2:size(inds,1)
    ref_roi=env(zinds,xinds,inds(i,1));
    cc=zeros(1,length(srch));
    for j=1:length(srch)
        tgt_roi=env(zinds,xinds+srch(j),inds(i,2));
        cc(j)=ncc(ref_roi,tgt_roi);
    end
    
    [lag,rho]=subsample_peak(cc,srch,'reconstructive_iterative');
    Y(i)=lag*dx;
    W(i)=rho; % can apply atanh(rho) here
end

% Weighted least squares
Wd=diag(real(W));
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
