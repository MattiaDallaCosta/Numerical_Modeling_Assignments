function PlotDisp(dispersion_ratio, h_over_lambda, nmodes, Data)
    % reshape: rows = modes, cols = different meshes
    % From a 1d matrix obtain a 2d matrix in which each column rapresents a
    % different mode
    R = reshape(dispersion_ratio, nmodes, []);
    H = reshape(h_over_lambda,  nmodes, []);

    figure; grid on; hold on;
    % Itarate for each line, which means each mode.
    for m = 1:nmodes
        [Hs, idx] = sort(H(m,:));
        plot(Hs, R(m,idx), '-o', 'LineWidth', 1.2, 'MarkerSize', 6);
    end
    %set(gca,'XScale','log')
    xlabel('h/\lambda');
    ylabel('k_h/k');
    title('Dispersion curves (per mode)');
    legend(arrayfun(@(m) sprintf('mode %d', m), 1:nmodes, 'UniformOutput', false), ...
           'Location','best');
    
    if (Data.save_sol_images)
        exportgraphics(gcf, fullfile('Plots', 'Dispersionplot.png'), 'Resolution', 300);
    end
end