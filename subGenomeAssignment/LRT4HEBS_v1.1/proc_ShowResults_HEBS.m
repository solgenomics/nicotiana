% Script to show a histogram and scatter plots of HEBS results

% Change default figure settings for increased aesthetics
set(0,'DefaultFigurePosition', [0 0 1200 600])
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultTextInterpreter','latex')
close all; clc;

% Bin width for histograms
BINWIDE = 0.25;

% Get all testable results and significant results
idx = isfinite(W);
idxBias = (W>Wadj);
idxUp = idxBias & HEBS_Results.HEBS > 0;
idxDn = idxBias & HEBS_Results.HEBS < 0;

% Determine number and mean HEBS for each group
Nup = sum(idxUp);
Ndn = sum(idxDn);
Mup = mean(HEBS_Results.HEBS(idxUp));
Mdn = mean(HEBS_Results.HEBS(idxDn));

% Make the plots
figure
subplot(1,2,1)
    histogram(HEBS_Results.HEBS(idx), 'binwidth', BINWIDE, 'facecolor', [0.4 0.4 0.4]); hold on
    histogram(HEBS_Results.HEBS(idxUp), 'binwidth', BINWIDE, 'facecolor', 'b');
    histogram(HEBS_Results.HEBS(idxDn), 'binwidth', BINWIDE, 'facecolor', 'y');
    xlabel('HEBS: Positive=Shift Towards Other in Leaf')
    ylabel('Count')
    xlim([-8 8])
    
    str = {'HEBS in M. luteus; Leaf vs. Petal', ... 
            ['Oth Bias: ' num2str(Nup) '; Mean HEBS: ' num2str(Mup)], ... 
            ['Gut Bias: ' num2str(Ndn) '; Mean HEBS: ' num2str(Mdn)]}';
    title(str)

subplot(1,2,2)
    scatter(F_Pars.HEB(idx), L_Pars.HEB(idx), '.k'); hold on
    scatter(F_Pars.HEB(idxUp), L_Pars.HEB(idxUp), 'ob', 'filled', 'MarkerEdgeColor', 'k');
    scatter(F_Pars.HEB(idxDn), L_Pars.HEB(idxDn), 'oy', 'filled', 'MarkerEdgeColor', 'k');
    xlabel('HEB in petal'); ylabel('HEB in leaf');
    legend('All data', 'HEBS significant and > 0', 'HEBS significant and < 0', 'location', 'northwest')

clearvars str idx* Nup Ndn Mup Mdn BINWIDE