% This example requires kWave
if ~exist('kWaveGrid', 'class')
  warning('kWave must be on the path to run this example.');
end

% Setup a system
grid = ScanCartesian('x', 1e-3*linspace(-20, 20, 1+40*2^3),'z', 1e-3*linspace(-02, 58, 1+60*2^3));
xdc = TransducerArray.L12_3v();
% xdc = TransducerArray('numel', 16, 'fc', 3e6, 'bw', [1.5, 4.5]*1e6, 'pitch', 1.5e3/3e6);
seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'c0', 1500);
us = UltrasoundSystem('scan', grid, 'xdc', xdc, 'seq', seq, 'fs', 10*xdc.fc);

kwargs = {
    'us', us;
    'grid', grid;
    'rho0', us.seq.c0 / 1.5;
    'c0', us.seq.c0;
    'scat_per_cell', 25;
    'ampdB', 0
};

% Convert kwargs to a structure for easier access
kwargsStruct = cell2struct(kwargs(:, 2), kwargs(:, 1), 1);

grid.dx = us.lambda / 4;
grid.dz = us.lambda / 4;
[Xg, ~, Zg] = grid.getImagingGrid();

% medium params
[rho0, c0] = deal(us.seq.c0 / 1.5, us.seq.c0);

% number of voxels per dim
wvl = range([grid.xb; grid.yb; grid.zb], 2) ./ us.lambda;
S = kwargsStruct.scat_per_cell * prod(wvl(wvl > 0));

% generate uniform random positions and normally distributed amplitudes
N   = arrayfun(@(p) grid.("n"+p), lower(grid.order)); % size in each dimension
ind = arrayfun(@(N) {randi(N, [S,1], 'like', rho0([]))}, N); % position indices
as  =                rand(    [S,1], 'like', rho0([]));      % amplitude values

% center point lesion assignment
zp = 1e-3*(10:10:50);
xp = 1e-3*(-15:7.5:15);
[Xp, Zp] = ndgrid(xp, zp);
centerPoints = [Xp(:), Zp(:)]; %  0*Xp(:) apply if y coordinate is applicable

% assign perturbation
rho = zeros(grid.size, 'like', rho0);
rho(sub2ind(grid.size, ind{:})) = as;


% add to base
rho = rho0 + db2pow(kwargsStruct.ampdB) .* rho0 * rho;
c = c0 * ones(grid.size);

% create a speckle perturbation distribution
drho = rho; % speckle perturbation distribution - not actually 0 

% add contrast lesions
add_lesion = true; % whether to add
ldB = -20; % lesion amplitude in dB
lesion_params = 1e-3 * [-12.5 10 3; 0 10 3; 12.5 10 3;...
    -12.5 20 3; 0 20 3; 12.5 20 3;...
    -12.5 30 3; 0 30 3; 12.5 30 3;...
    -12.5 40 3; 0 40 3; 12.5 40 3;...
    -12.5 50 3; 0 50 3; 12.5 50 3]; % (x0, z0, r) mm e.g. lateral, axial, radius
if add_lesion
    for l = 1:size(lesion_params, 1)
        % find the indices within the lesion
        p = lesion_params(l,:); % (x0, z0, r)

        % logical Mask for points within the lesion
        i = (Xg - p(1)) .^ 2 + (Zg - p(2)).^2 < p(3)^2;
        
        % modify the density perturbation distribution within the lesion
        drho(i) = drho(i) .* db2pow(ldB); % #ok<SAGROW> % not actually 0
    end
end

% add the perturbation
rho = rho + drho;

% sample the Medium on the grid
% args = namedargs2cell(rmfield(kwargsStruct, ["scat_per_cell","ampdB"]));
med = Medium.Sampled(grid, c, rho, 'c0',c0, 'rho0',rho0);

figure; 
imagesc(med, us.scan, 'props', ["c","rho"]);

% Simulate the ChannelData
us.fs = single(us.fs); % accelerate
chd = kspaceFirstOrder(us, med, grid, 'CFL_max', 0.3); % CLF_max lower val to 0.25
% sct = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
% chd = greens(us, sct);

% cast to gpu array to circumnavigate ram clutter
chd = gpuArray(chd);

% convert chd to complex
chd = hilbert(chd);

% run beamformers
b_naive = DAS(us,chd);
b_c0    = bfEikonal(us, chd, med, grid);
b_axial = bfAxialGlobalAvg(us,chd,med,grid);

figure;
imagesc(hilbert(chd));
title('Channel Data');

figure;
bim_c0    = mod2db(b_c0   );
imagesc(us.scan, bim_c0   , [-40 0] + max(bim_c0(:)   ));
colormap gray; colorbar;
title('Eikonal Delay-and-Sum Contrast Lesions');

figure;
bim_c0    = mod2db(b_axial   );
imagesc(us.scan, bim_c0   , [-40 0] + max(bim_c0(:)   ));
colormap gray; colorbar;
title('AGA Delay-and-Sum Contrast Lesions');

figure;
bim_c0    = mod2db(b_naive   );
imagesc(us.scan, bim_c0   , [-40 0] + max(bim_c0(:)   ));
colormap gray; colorbar;
title('Naive Delay-and-Sum Contrast Lesions');

% Save workspace
save('contrast_lesion_medium');
