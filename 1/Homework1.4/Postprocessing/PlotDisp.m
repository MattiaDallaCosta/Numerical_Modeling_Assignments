function PlotDisp(dispersion_ratio, h_over_lambda, nmodes, Data)
    % reshape: rows = modes, cols = different meshes
    % From a 1d matrix obtain a 2d matrix in which each column rapresents a
    % different mode
    R = reshape(dispersion_ratio, nmodes, []);
    H = reshape(h_over_lambda,  nmodes, []);

    fs        = 30;
    lw        = 1.9;
    fnt       = 'Times New Roman';
    interp_tex = 'latex';
    
    clr = [ ...
        0.00  0.45  0.74;   %  1  blue
        0.85  0.33  0.10;   %  2  orange
        0.47  0.67  0.19;   %  3  green
        0.49  0.18  0.56;   %  4  purple
        0.93  0.69  0.13;   %  5  yellow
        0.30  0.75  0.93;   %  6  cyan
        0.64  0.08  0.18;   %  7  dark red
        0.25  0.25  0.25;   %  8  dark gray
        0.00  0.60  0.50;   %  9  teal
        0.75  0.00  0.75;   % 10  magenta
        0.60  0.60  0.60];  % 11  light gray

    fig4 = figure('Units','centimeters','Position',[2 2 30 16], ...
              'Color','w','PaperPositionMode','auto');
    hold on; box on; grid on;
    set(gca, 'GridLineStyle','-', 'GridAlpha',0.15);
    
    for m = 1:nmodes
        [Hs, idx] = sort(H(m,:));
        ci = clr(mod(m-1, size(clr,1)) + 1, :);
        plot(Hs, R(m,idx), '-o', ...
             'LineWidth',lw, 'MarkerSize',6, 'Color',ci);
    end
    
    set(gca, 'XScale','log', 'YScale','log');
    
    xlabel('$h/\lambda$', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    ylabel('$k_h/k$',     'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    title('Dispersion curve', 'Interpreter',interp_tex, 'FontName',fnt, 'FontSize',fs);
    
    legend(arrayfun(@(m) sprintf('$\\mathrm{mode}\\ %d$', m), ...
                    1:nmodes, 'UniformOutput', false), ...
           'Interpreter',interp_tex, ...
           'Location','best', ...
           'FontName',fnt, ...
           'FontSize',fs-8, ...
           'Box','on');
    
    set(gca, 'FontName',fnt, 'FontSize',fs-4, 'TickLabelInterpreter','latex');
    
    if Data.save_sol_images
        exportgraphics(fig4, fullfile('Plots', 'Dispersionplot.pdf'), 'ContentType','vector');
    end
end