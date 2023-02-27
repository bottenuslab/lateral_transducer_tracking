% [lag,rho] = subsample_peak(cc, lags, method, (parameters))
% 
% Perform subsample peak estimation of a cross-correlation curve using
% varoius algorithms.
% cc - Cross-correlation curve
% lags - Corresponding lags for cross-correlation curve (integer spacing)
% method - 'spline','cubic','parabolic','cosine','reconstructive','reconstructive_iterative'
%
% See the code for individual methods for specific paremeters passed in as
% name/value pairs
%
% Reconstructive interpolation implemented from Cespedes et al, Ultrasonic
% Imaging, 1995
%
function [lag,rho] = subsample_peak(cc,lags,method,varargin)

p=inputParser;
cc=cc(:)'; % Make row vectors
lags=lags(:)';

switch method
    
    case 'spline'
        addParameter(p,'upsample',100); % Factor to upsample estimate by
        p.parse(varargin{:});
        upsample=p.Results.upsample;
        
        lags_up=lags(1):1/upsample:lags(end);
        cur=interp1(lags,cc,lags_up,'spline');
        lag=mean(lags_up(cur==max(cur)));
        rho=max(cur);
        
    case 'cubic'
        addParameter(p,'upsample',100); % Factor to upsample estimate by
        p.parse(varargin{:});
        upsample=p.Results.upsample;
        
        lags_up=lags(1):1/upsample:lags(end);
        cur=interp1(lags,cc,lags_up,'cubic');
        lag=mean(lags_up(cur==max(cur)));
        rho=max(cur);
        
    case 'parabolic'
        ind=find(cc==max(cc),1,'first');
        if(ind==1||ind==length(cc)) % Don't subsample if at the end
            lag=lags(ind);
            rho=max(cc);
        else
            cur=cc(ind-1:ind+1);
            offset=1/2*(cur(1)-cur(3))/(cur(1)-2*cur(2)+cur(3));
            lag=lags(ind)+offset;
            rho=(cur(1)-2*cur(2)+cur(3))/2*offset.^2+(cur(3)-cur(1))/2*offset+cur(2);
        end
        
    case 'cosine'
        ind=find(cc==max(cc),1,'first');
        if(ind==1||ind==length(cc)) % Don't subsample if at the end
            lag=lags(ind);
            rho=max(cc);
        else
            cur=cc(ind-1:ind+1);
            w0=acos((cur(1)+cur(3))/(2*cur(2)));
            theta=atan((cur(1)-cur(3))/(2*cur(2)*sin(w0)));
            offset=-theta/w0;
            lag=lags(ind)+offset;
            rho=0; %??
        end
        
    case 'reconstructive'
        addParameter(p,'len',4); % Half-width of window to fit
        addParameter(p,'upsample',100); % Factor to upsample estimates by
        p.parse(varargin{:});
        len=p.Results.len;
        upsample=p.Results.upsample;
        
        % Find peak
        ind=find(cc==max(cc));
        % Find offset around peak
        len=min([len,ind-1,length(cc)-ind]);
        search=-len:len;
        search_up=search(1):1/upsample:search(end);
        cc_up=sum(cc(ind+search).*sinc(search_up'-search),2);
        lag=lags(ind)+search_up(cc_up==max(cc_up));
        rho=max(cc_up);
        
    case 'reconstructive_iterative'
        addParameter(p,'iter',10); % Number of iterations
        addParameter(p,'len',4); % Half-width of window to fit
        p.parse(varargin{:});
        iter=p.Results.iter;
        len=p.Results.len;
        
        % Find peak
        ind=find(cc==max(cc),1,'first');
        if(ind==1||ind==length(cc))  % Don't subsample estimate if point is at the end
            lag=lags(ind);
            rho=max(cc);
        else
            % Set interpolation region around peak
            len=min([len,ind-1,length(cc)-ind]);
            search=-len:len;
            % Find largest 3 values
            cur=cc(ind-1:ind+1);
            if(cur(1)>cur(3)) % work from left
                locs=[0 1];
                for j=1:iter
                    % Interpolate middle value
                    mid=mean(locs);
                    cc_up=sum(cc(ind+search).*sinc(mid-search));
                    if(cc_up>cur(1)) % Step in the correct direction
                        locs=[mid locs(2)];
                    else
                        locs=[locs(1) mid];
                    end
                end
                lag=lags(ind)+mean([-1 mid]); % Average left and right points
            else % work from right
                locs=[-1 0];
                for j=1:iter
                    % Interpolate middle value
                    mid=mean(locs);
                    cc_up=sum(cc(ind+search).*sinc(mid-search));
                    if(cc_up>cur(3)) % Step in the correct direction
                        locs=[locs(1) mid];
                    else
                        locs=[mid locs(2)];
                    end
                end
                lag=lags(ind)+mean([mid 1]); % Average left and right points
            end
            rho=sum(cc(ind+search).*sinc(lag-lags(ind)-search)); % Get correlation at peak
        end
    
end


