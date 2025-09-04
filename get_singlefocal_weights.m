function [wnew, psuedo_deg] = get_singlefocal_weights(beamAng,beamdist, N, fc, whichbeam, wannaplot)
% Ish Jain
% Dec 15 2020
% 
% EDIT: Luke Wilson Sept 2025
%       I modified Ish's get_null_singlebeam_weights to work as a beam
%       focuser weight generator instead. given the angle and distance of
%       the focal point, it will point a beam in the proper direction and
%       gaurantee phase alignment at the focal point
%
%
%

c = 3e8;
lambda = c/fc;
d_ant = lambda/2;


if(nargin<6)
    wannaplot=1;
end
if(nargin<5)
    whichbeam=1; %whichbeam in beamAng is main beam, others are null
end
if(nargin<1)
    beamAng = [0,10];
end

% beam directions and distance
beamAOD = beamAng(whichbeam); % main beam direction
nullidx=setdiff(1:length(beamAng),whichbeam);
beamNull = beamAng(nullidx); %null beam direction

% first do phase allignment for focal distance
poi_x = beamdist * sind(beamAOD);
poi_y = beamdist * cosd(beamAOD);

distances_for_phases = zeros(1,N); 
phases_for_weights = zeros(1,N);

    ant_x = 0:d_ant:(N-1)*d_ant;
    ant_x = ant_x - mean(ant_x);
    ant_y = zeros(1,N);

for i=1:N

    dx = poi_x - ant_x(i);
    dy = poi_y - ant_y(i);

    distances_for_phases(i) = sqrt(dx.^2 + dy.^2);
    phases_for_weights(i) = exp(1j * 2*pi .*  ...
        distances_for_phases(i) ./ lambda);
end

% At this point I have the phases to get allignment,
% now i need to see what the angle effect is by treating it
% as a "psuedo-weight-vector" and solving for psi

% General kd with arbitrary spacing d_ant (already λ/2 in your code)
k  = 2*pi/lambda;
kd = k * d_ant;

psi_steps = unwrap(angle(phases_for_weights(2:end) .* conj(phases_for_weights(1:end-1))));
psi_hat   = mean(psi_steps);              % average linear phase step
sin_th    = psi_hat / kd;                 % kd * sinθ ≈ ψ
sin_th    = max(min(sin_th,1),-1);
psuedo_deg = asind(sin_th);

% Compensate the requested beam angle by removing the linear component
% beamAng = beamAng - psuedo_deg;    % works for positive/negative naturally


%Number of antennas in ULA
theta = [-90:0.1:90].'; %Range of theta for plotting beampattern

for bid=1:length(beamAng) %create standard single beam patterns
    w(:,bid) = get_multibeam_weights(beamAng(bid),1,0,N);
end

%create nulling
psiall = 2*pi/2*sind(theta); % grid of psi for plotting
for n=1:N
    V(n,:) =  exp(1j*(n-1)*psiall); % array manifold vector
end


%% beam pattern after applying null constraint

%Apply null constraint and get weights
vc = w(:,nullidx);%%V(:,psindx); 
%  vc=vc.';
% Pc = vc*vc'/(vc'*vc); %vector computation
Pc = vc/(vc'*vc)*vc'; %matrix computation: I did this computation on my notebook
Pcp = eye(N)-Pc;
wnew= Pcp*w(:,whichbeam);

% Return beam information
bm.Bsingle = V.'*wnew;
bm.theta=theta;

% Normalize
wnew = wnew./norm(wnew);
if(wannaplot)
    figure(111); clf;
    plot_beam_pattern(wnew, bm.theta,N)
end
end

