function avg_error = target_points_eval(rho,bimage)

    %%%%%
    % get true point target location
    %%%%%
    [row,col]   = find(rho == 2000); % modify density scatter val
    ptIndices   = [row,col];

    %%%%%
    % Thresholding
    %%%%%
    flattened_bim   = bimage(:);
    thr_bim         = max(flattened_bim) * 0.75; % modify threshold val
    [rows, cols]    = size(bimage);
        
    thresh_bim = bimage;
    for r = 1:rows
        for c = 1:cols
            if thresh_bim(r,c) <= thr_bim
                thresh_bim(r,c) = 60;   % modify dB for dismissed points
            end
         end
    end
    

    %%%%%
    % error calc using euclidean distance
    %%%%%
    error_sum = 0;
    for pt = 1:length(ptIndices)
        sim_ptIndx = concentric_square_search(thresh_bim, [ptIndices(pt,1),ptIndices(pt,2)]);
        error_sum = error_sum + norm(sim_ptIndx-[ptIndices(pt,1),ptIndices(pt,2)]);
    end
    disp(error_sum * (1/length(ptIndices)));
    avg_error = error_sum * (1/length(ptIndices));
end
