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
        plot(Hs, R(m,idx), '-o', 'LineWidth', 1.7, 'MarkerSize', 6);
    end
    %set(gca,'YScale','log')
    %set(gca,'XScale','log')
    xlabel('h/\lambda', fontsize=17);
    ylabel('k_h/k', fontsize=17);
    title('Dispersion curves (per mode)', fontsize=18);
    legend(arrayfun(@(m) sprintf('mode %d', m), 1:nmodes, 'UniformOutput', false), ...
           'Location','best', fontsize=18);
    
    if (Data.save_sol_images)
        exportgraphics(gcf, fullfile('Plots', 'Dispersionplot.png'), 'Resolution', 300);
    end
end