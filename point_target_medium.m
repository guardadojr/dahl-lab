% This example requires kWave
if ~exist('kWaveGrid', 'class')
  warning('kWave must be on the path to run this example.');
end

% Setup a system
sscan = ScanCartesian('x', 1e-3*linspace(-20, 20, 1+40*2^3),'z', 1e-3*linspace(-02, 58, 1+60*2^3));
xdc = TransducerArray.L12_3v();
% xdc = TransducerArray('numel', 16, 'fc', 3e6, 'bw', [1.5, 4.5]*1e6, 'pitch', 1.5e3/3e6);
seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'c0', 1500);
us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'seq', seq, 'fs', 10*xdc.fc);
% [us.scan.dx, us.scan.dz] = deal(us.lambda / 2);

%
zp = 1e-3*(10:10:50);
xp = 1e-3*(-15:7.5:15);
[Xp, Zp] = ndgrid(xp, zp);
pp = [Xp(:), 0*Xp(:), Zp(:)];

sct = Scatterers('pos', pp');
sct.c0 = us.seq.c0;

% Create a Medium to simulate
[c0, rho0] = deal(1.5e3, 1e3);
[c, rho] = deal(c0*ones(sscan.size), rho0*ones(sscan.size));
[Xg, ~, Zg] = sscan.getImagingGrid();

% Define isoimpedance layers
z0 = rho0 * c0; % ambient impedance
[c(Zg > 15e-3), rho(Zg > 15e-3)] = deal(1.4e3, z0/1.4e3); % isoimpedance
[c(Zg > 25e-3), rho(Zg > 25e-3)] = deal(1.6e3, z0/1.6e3); % isoimpedance
[c(Zg > 35e-3), rho(Zg > 35e-3)] = deal(1.4e3, z0/1.4e3); % isoimpedance
[c(Zg > 45e-3), rho(Zg > 45e-3)] = deal(1.5e3, z0/1.5e3); % isoimpedance

%{ 
% Define density scatterers
rho(Xg == 0 & Zg == 10e-3) = rho0*2;
rho(Xg == 0 & Zg == 20e-3) = rho0*2; 
rho(Xg == 0 & Zg == 30e-3) = rho0*2; 
rho(Xg == 0 & Zg == 40e-3) = rho0*2; 
rho(Xg == 0 & Zg == 50e-3) = rho0*2;

rho(Xg == -12.5e-3 & Zg == 10e-3) = rho0*2;
rho(Xg == -12.5e-3 & Zg == 20e-3) = rho0*2; 
rho(Xg == -12.5e-3 & Zg == 30e-3) = rho0*2; 
rho(Xg == -12.5e-3 & Zg == 40e-3) = rho0*2; 
rho(Xg == -12.5e-3 & Zg == 50e-3) = rho0*2;

rho(Xg == 12.5e-3 & Zg == 10e-3) = rho0*2;
rho(Xg == 12.5e-3 & Zg == 20e-3) = rho0*2; 
rho(Xg == 12.5e-3 & Zg == 30e-3) = rho0*2; 
rho(Xg == 12.5e-3 & Zg == 40e-3) = rho0*2; 
rho(Xg == 12.5e-3 & Zg == 50e-3) = rho0*2;
%}

% add the density scatterers
pg = sscan.positions(); % 3 x Z x X x Y
for s = 1:sct.numScat
    % find the closest grid point
    [~, i] = min(vecnorm(pg - sct.pos(:,s),2,1),[],'all','linear');

    % double the density at that grid point
    rho(i) = rho(i) * 2;
end

%{
% create a speckle perturbation distribution
drho = 0; % speckle perturbation distribution - not actually 0 

% add contrast lesions
add_lesion = true; % whether to add
ldB = -20; % lesion amplitude in dB
lesion_params = 1e-3 * [0 30 5; 10 40 5]; % (x0, z0, r) mm e.g. lateral, axial, radius
if add_lesion
    for l = 1:size(lesion_params, 1)
        % find the indices within the lesion
        p = lesion_params(l,:); % (x0, z0, r)
        i = (Xg - p(1)) .^ 2 + (Zg - p(2)).^2 < p(3);
        
        % modify the density perturbation distribution within the lesion
        drho(i) = drho(i) .* db2pow(ldB); %#ok<SAGROW> % not actually 0
    end
end

% add the perturbation
rho = rho + drho;

%}

% Construct the Medium
med = Medium.Sampled(sscan, c, rho, 'c0', c0, 'rho0', rho0);

figure; 
imagesc(med, us.scan, 'props', ["c","rho"]);

%% Simulate the ChannelData
us.fs = single(us.fs); % accelerate
chd = kspaceFirstOrder(us, med, sscan, 'CFL_max', 0.5);
% chd = greens(us, sct);

% Beamform
b_naive = DAS(us, chd, 'apod', 1);
b_c0    = bfEikonal(us, chd, med, sscan, 'apod', 1);
b_axial = bfAxialGlobalAvg(us,chd,med,sscan);

% Display the ChannelData
figure;
imagesc(hilbert(chd));
title('Channel Data');

% Display the images
bim_naive = mod2db(b_naive);
bim_c0    = mod2db(b_c0   );
bim_axial = mod2db(b_axial);

% Create a figure with two subplots and save as one PNG
figure;
imagesc(us.scan, bim_naive, [-40 0] + max(bim_naive(:)));
colormap gray; colorbar;
title('Naive Delay-and-Sum');

figure;
imagesc(us.scan, bim_c0   , [-40 0] + max(bim_c0(:)   ));
colormap gray; colorbar;
title('Eikonal Delay-and-Sum');

figure;
imagesc(us.scan, bim_axial   , [-40 0] + max(bim_axial(:)   ));
colormap gray; colorbar;
title('AGA Delay-and-Sum');

% Save workspace
save('point_target_medium');
