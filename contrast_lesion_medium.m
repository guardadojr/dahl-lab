% This example requires kWave
if ~exist('kWaveGrid', 'class')
  warning('kWave must be on the path to run this example.');
end

% Setup a system
sscan = ScanCartesian('x', 1e-3*linspace(-20, 20, 1+40*2^3),'z', 1e-3*linspace(-02, 58, 1+60*2^3));
% xdc = TransducerArray.L12_3v();
xdc = TransducerArray('numel', 16, 'fc', 3e6, 'bw', [1.5, 4.5]*1e6, 'pitch', 1.5e3/3e6);
seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'c0', 1500);
us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'seq', seq, 'fs', 10*xdc.fc);

% medium params
wvl = range([grid.xb; grid.yb; grid.zb], 2) ./ us.lambda;
S = kwargs.scat_per_cell * prod(wvl(wvl > 0));

% generate uniform random positions and normally distributed amplitudes
N   = arrayfun(@(p) grid.("n"+p), lower(grid.order)); % size in each dimension
ind = arrayfun(@(N) {randi(N, [S,1], 'like', rho0([]))}, N); % position indices
as  =                rand(    [S,1], 'like', rho0([]));      % amplitude values

% assign perturbation
rho = zeros(grid.size, 'like', rho0);
rho(sub2ind(grid.size, ind{:})) = as;

% add to base
rho = rho0 + db2pow(kwargs.ampdB) .* rho0 * rho;
c = c0 * ones(grid.size);

% sample the Medium on the grid
args = namedargs2cell(rmfield(kwargs, ["scat_per_cell","ampdB"]));
med = Medium.Sampled(grid, c, rho, args{:});

