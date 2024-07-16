
% This example requires kWave
  if ~exist('kWaveGrid', 'class')
      warning('kWave must be on the path to run this example.');
  end
  
  % Setup a system
  sscan = ScanCartesian('x', 1e-3*linspace(-20, 20, 1+40*2^3),'z', 1e-3*linspace(-02, 58, 1+60*2^3));
  xdc = TransducerArray('numel', 16, 'fc', 3e6, 'bw', [1.5, 4.5]*1e6, 'pitch', 1.5e3/3e6);
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
  
  % Define density scatterers
  rho(Xg == 0 & Zg == 10e-3) = rho0*2;
  rho(Xg == 0 & Zg == 20e-3) = rho0*2; 
  rho(Xg == 0 & Zg == 30e-3) = rho0*2; 
  rho(Xg == 0 & Zg == 40e-3) = rho0*2; 
  rho(Xg == 0 & Zg == 50e-3) = rho0*2;
  
  % Construct the Medium
  med = Medium.Sampled(sscan, c, rho);
  
  % Simulate the ChannelData
  us.fs = single(us.fs); % accelerate
  chd = kspaceFirstOrder(us, med, sscan, 'CFL_max', 0.5);
  
  % Beamform
  b_naive = DAS(us, chd);
  b_c0    = bfEikonal(us, chd, med, sscan);
  
  % Display the ChannelData
  figure;
  imagesc(real(chd));
  title('Channel Data')
  
  % Display the images
  bim_naive = mod2db(b_naive);
  bim_c0    = mod2db(b_c0   );
  
  figure;
  subplot(1,2,1);
  imagesc(us.scan, bim_naive, [-80 0] + max(bim_naive(:)));
  colormap gray; colorbar;
  title('Naive Delay-and-Sum');
  
  subplot(1,2,2);
  imagesc(us.scan, bim_c0   , [-80 0] + max(bim_c0(:)   ));
  colormap gray; colorbar;
  title('Eikonal Delay-and-Sum');





  %us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'seq', seq, 'fs', 10*xdc.fc);
  %loactionTP1 = concentric_square_search(thresh_bim_c0, [334,161]);
  %target_points_eval(rho,bim_naive)

 %{
  %%%%%%%%%%%
  % increase resolution
  %%%%%%%%%%%
  [sscan.dx, sscan.dz] = deal(us.lambda/12);

  % Beamform
  b_naive_HR = DAS(us, chd);
  b_c0_HR    = bfEikonal(us, chd, med, sscan);
  
  % Display the ChannelData
  figure;
  imagesc(real(chd));
  title('Channel Data Higher Resolution')
  
  % Display the images
  bim_naive_HR = mod2db(b_naive_HR);
  bim_c0_HR    = mod2db(b_c0_HR   );
  
  figure;
  subplot(1,2,1);
  imagesc(us.scan, bim_naive_HR, [-80 0] + max(bim_naive_HR(:)));
  colormap gray; colorbar;
  title('Naive Delay-and-Sum Higher Resolution');
  
  subplot(1,2,2);
  imagesc(us.scan, bim_c0_HR   , [-80 0] + max(bim_c0_HR(:)   ));
  colormap gray; colorbar;
  title('Eikonal Delay-and-Sum Higher Resolution');
  %}
