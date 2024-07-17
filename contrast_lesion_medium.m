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