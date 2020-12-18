%% Add paths

%--------------------
% GrOpt Matlab Libraries
% Download with: git clone https://github.com/cmr-group/gropt
%--------------------
addpath(genpath('../gropt/matlab'));

%--------------------
% Pulseq Libraries
% Download with: git clone https://github.com/pulseq/pulseq
%--------------------
addpath(genpath('../pulseq/matlab'));


%% Start Sequence 
% Code is based off of the Pulseq demo "writeEpiSpinEcho.m"

seq=mr.Sequence();              % Create a new sequence object
fov=256e-3; Nx=64; Ny=64;       % Define FOV and resolution
TE=80e-3;

gmax = 32; % Gradient Amplitude limit [mT/m]
smax = 120; % Slew Rate Limit [T/m/s]

% Set system limits for GrOpt
gropt_params = struct;
gropt_params.mode = 'diff_beta'; % Add beta maximization to GrOpt for diffusion
gropt_params.gmax = 0.99*gmax;
gropt_params.smax = 0.99*smax;
gropt_params.dt = 40e-6; % Raster time to optimize at (>10us for speed)
gropt_params.dt_out = 10e-6;  % Gradient raster time

% Set system limits for Pulseq
lims = mr.opts('MaxGrad', gmax, 'GradUnit', 'mT/m',...
    'MaxSlew', smax, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  

% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,lims,'Duration',3e-3,...
    'SliceThickness',3e-3,'apodization',0.5,'timeBwProduct',4);

% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
readoutTime = 3.2e-4;
gx = mr.makeTrapezoid('x',lims,'FlatArea',kWidth,'FlatTime',readoutTime);
adc = mr.makeAdc(Nx,lims,'Duration',gx.flatTime,'Delay',gx.riseTime);

% Pre-phasing gradients
preTime=8e-4;
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);

% Readout prewinders
RO_preTime=4e-4;
gx_readout_pre = mr.makeTrapezoid('x',lims,'Area',-(gx.area/2-deltak/2),'Duration',RO_preTime);
gy_readout_pre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak,'Duration',RO_preTime);

% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur);

% Refocusing pulse with spoiling gradients
rf180 = mr.makeBlockPulse(pi,lims,'Duration',500e-6,'use','refocusing');
[rf180, gz180] = mr.makeSincPulse(pi,lims,'Duration',3e-3,...
                'SliceThickness',8e-3,'apodization',0.5,'timeBwProduct',4,'use','refocusing');
gzSpoil = mr.makeTrapezoid('z',lims,'Area',gz.area*2,'Duration',3*preTime);


% Calculate delay time
durationToCenter = (Nx/2+0.5)*mr.calcDuration(gx) + Ny/2*mr.calcDuration(gy);
rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);
rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);

% ------------------------
% GrOpt Waveforms
% ------------------------
% Calculate 90 and 180 deadtimes
T_90 = gz.fallTime + gz.flatTime/2 + preTime;
T_180 = 6*preTime + 3e-3;

% Time needed for prewinders and EPI to reach TE (more dead time for diffusion) 
T_readout = durationToCenter + RO_preTime;

gropt_params.MMT = 0;  % Null the 0th moments (area)
gropt_params.T_readout = T_readout;   % Set timings 
gropt_params.T_90 = T_90;  
gropt_params.T_180 = T_180; 
gropt_params.TE = TE;
gropt_params.eddy_params = 60; % Add constraint to null eddy current time constant 60ms

% Run GrOpt to get the waveform
G_diff = gropt(gropt_params);

% Split waveform at 180 and remove zeros so we can use as 2 individual
% blocks
[G1, G2] = split_diff(G_diff, gropt_params);

% Make Pulseq arbitrary gradients from the GrOpt split waveforms
gdiff1 = mr.makeArbitraryGrad('x',lims.gamma*G1, 'system',lims);
gdiff2 = mr.makeArbitraryGrad('x',lims.gamma*G2, 'system',lims);


% Define sequence blocks
seq.addBlock(rf,gz);
seq.addBlock(gzReph); 
seq.addBlock(gdiff1);
seq.addBlock(gzSpoil);
seq.addBlock(rf180, gz180);
seq.addBlock(gzSpoil);
seq.addBlock(gdiff2);
seq.addBlock(gx_readout_pre, gy_readout_pre); 
for i=1:Ny
    seq.addBlock(gx,adc);           % Read one line of k-space
    seq.addBlock(gy);               % Phase blip
    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
end
seq.addBlock(mr.makeDelay(1e-4));

seq.write('gropt_encode.seq');   % Output sequence for scanner
seq.plot();                % Plot sequence waveforms


%% Plot diffusion waveform and Eddy Current spectrum

% Plot gradient waveform
tt = [0:numel(G_diff)-1] .* gropt_params.dt_out .* 1e3;
figure();
set(gcf,'Position',[100 100 900 500])
plot(tt, G_diff.*1000, 'LineWidth', 2);
hold on;
plot([min(xlim()),max(xlim())],[0,0], 'k--');
title('GrOpt Waveform'); xlabel('t [ms]'); ylabel('G [mT/m]');
ax = gca;
ax.FontSize = 18; 


% Plot eddy currents
lam_max = 120; % maximum lambda [ms] for eddy spectrum analysis
[all_lam, all_e] = get_eddy(lam_max, G_diff, gropt_params.dt_out);
figure();
set(gcf,'Position',[100 100 900 500])
plot(all_lam, all_e, 'LineWidth', 2, 'Color', '#d81b60');
hold on;
plot([min(xlim()),max(xlim())],[0,0], 'k--');
title('Eddy Current Spectrum'); xlabel('lambda [ms]');
ax = gca;
ax.FontSize = 18; 
