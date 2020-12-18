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

%% Make Sequence
% Code is based off of the Pulseq GRE seqeunce


seq=mr.Sequence();              % Create a new sequence object
fov=360e-3; Nx=128; Ny=128;      % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=8e-3;            % slice


gmax = 32; % Gradient Amplitude limit [mT/m]
smax = 150; % Slew Rate Limit [T/m/s]

% Set system limits for GrOpt
gropt_params = struct;
gropt_params.mode = 'free'; % Tell GrOpt this is a feasibility problem
gropt_params.gmax = 0.99*gmax;
gropt_params.smax = 0.99*smax;
gropt_params.dt = 10e-6;  % Gradient raster time


% Set system limits
lims = mr.opts('MaxGrad', gmax, 'GradUnit', 'mT/m', ...
    'MaxSlew', smax, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',2e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,...
    'system',lims);

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
deltak=1/fov;

% Readout
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',3.6e-3,'system',lims);
adc = mr.makeAdc(2*Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',lims);

% ------------------------
% GrOpt Waveforms
% ------------------------
gamma = 42.577478518;

% Get 1/2 of area of slice select to refocus 
M0_refocus = gz.area /2 / 42.577478518;

% Calculate M1 parameters for velocity encoding
venc = 0.8;  % Venc of 0.8 m/s 
delta_M1 = 1e3 / venc / 2 / gamma;  % M1 required for the Venc
shift_M1 = -3.8;  % M1 shifting for optimal TE

% Set GrOpt parameters
gropt_params.T = 1.3;  % Duration of bipolar
gropt_params.pns_thresh = 0.8; % Constrain the waveform to PNS < 0.8

% Design one bipolar ("flow up")
gropt_params.M0 = -M0_refocus;
gropt_params.M1 = delta_M1/2 + shift_M1;

% Get flow up waveform from GrOpt
G_up = gropt(gropt_params);


% Design second bipolar ("flow down")
gropt_params.M0 = -M0_refocus;
gropt_params.M1 = -delta_M1/2 + shift_M1;

% Get flow down waveform from GrOpt
G_down = gropt(gropt_params);


% Create Pulseq gradients from GrOpt Waveforms
g_flowup = mr.makeArbitraryGrad('z',lims.gamma*G_up, 'system',lims);
g_flowdown = mr.makeArbitraryGrad('z',lims.gamma*G_down, 'system',lims);


% To minimize TE we find the shortest time in x, y, and z for
% rephasing/prewinding
gx_mid0 = mr.makeTrapezoid('x','Area',-gx.area/2,'system',lims);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
gy_mid0 = mr.makeTrapezoid('z','Area',phaseAreas(1),'system',lims);

mid_dur = mr.calcDuration(gx_mid0, gy_mid0, g_flowup, g_flowdown);

% To minimize TR we also find the shortest time in x, y, and z for
% spoiling
gx_end0 = mr.makeTrapezoid('x','Area',2*Nx*deltak,'system',lims);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
gy_end0 = mr.makeTrapezoid('z','Area',phaseAreas(1),'system',lims);
gz_end0 = mr.makeTrapezoid('z','Area',4/sliceThickness,'system',lims);

end_dur = mr.calcDuration(gx_end0, gy_end0, gz_end0);


% Now make the gradients with those calculated times
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration', mid_dur, 'system',lims);
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak, 'Duration', end_dur, 'system',lims);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'Duration', end_dur,'system',lims);


% Put sequence together
rf_phase=0;
rf_inc=0;

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    for i_vel = 1:2 % velocity encodes 
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
   
        seq.addBlock(rf,gz);
        gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),...
                                 'Duration',mid_dur, 'system',lims);
        if i_vel == 1
            seq.addBlock(gxPre,gyPre,g_flowup);
        elseif i_vel == 2
            seq.addBlock(gxPre,gyPre,g_flowdown);
        end
        seq.addBlock(gx,adc);
        gySpoil = mr.makeTrapezoid('y','Area',-phaseAreas(i),...
                                 'Duration',end_dur, 'system',lims);
        seq.addBlock(gxSpoil,gySpoil,gzSpoil)
    end
end

%% Pulseq Plot
seq.write('gropt_pns_pcmri.seq');   % Output sequence for scanner
seq.plot();                      % Plot sequence waveforms

%% Plot waveforms and PNS prediction

% Plot gradient waveform
tt = [0:numel(G_up)-1] .* gropt_params.dt .* 1e3;
figure();
set(gcf,'Position',[100 100 900 500])
plot(tt, G_up.*1000, 'LineWidth', 2);
hold on;
plot(tt, G_down.*1000, 'LineWidth', 2);
plot([min(xlim()),max(xlim())],[0,0], 'k--');
title('GrOpt Waveform'); xlabel('t [ms]'); ylabel('G [mT/m]');
ax = gca;
ax.FontSize = 18; 


figure()
set(gcf,'Position',[100 100 900 500])
y = abs(get_stim(G_up, gropt_params.dt)');
plot(tt, y, 'LineWidth', 2, 'Color', '#00695c');
ylabel('PNS [A.U.]');
xlabel('t [ms]');
hline = refline([0 0.8]);
hline.Color = 'r';
hline.LineStyle = '--';
hline.LineWidth = 2;
title('Peripheral Nerve Stimulation');
ylim([0, 1.0])
ax = gca;
ax.FontSize = 18; 


