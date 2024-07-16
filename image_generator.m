function image_generator(bim_c0, bim_naive, us)
    figure;
    subplot(2,2,1);
    imagesc(us.scan, bim_naive, [-80 0] + max(bim_naive(:)));
    colormap gray; colorbar;
    title('Naive Delay-and-Sum');
    
    subplot(2,2,2);
    imagesc(us.scan, bim_c0   , [-80 0] + max(bim_c0(:)   ));
    colormap gray; colorbar;
    title('Eikonal Delay-and-Sum');

    thr_bim_naive = threshold(bim_naive);
    subplot(2,2,3);
    imagesc(us.scan, thr_bim_naive   , [-80 0] + max(thr_bim_naive(:)   ));
    colormap gray; colorbar;
    title('Thresh Naive Delay-and-Sum');

    thr_bim_c0 = threshold(bim_c0); 
    subplot(2,2,4);
    imagesc(us.scan, thr_bim_c0   , [-80 0] + max(thr_bim_c0(:)   ));
    colormap gray; colorbar;
    title('Thresh Eikonal Delay-and-Sum');

    function thresh_bim = threshold(bimage)
        %%%%%
        % Thresholding
        %%%%%
        flattened_bim   = bimage(:);
        thr_bim         = max(flattened_bim) * 0.75; % modify threshold val
        [rows, cols]    = size(bimage);

        for r = 1:rows
            for c = 1:cols
                if bimage(r,c) <= thr_bim
                    bimage(r,c) = 60;   % modify dB for dismissed points
                end
             end
        end
        thresh_bim = bimage;

    end
    
end