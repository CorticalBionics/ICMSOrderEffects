% ISIAmpDisc Analysis
addpath(genpath("HelperFunctions"))
load(fullfile(DataPath(), 'ISIAmpDiscData_Processed.mat'))
load(fullfile(DataPath(), 'MechAmpDiscrimData_Processed.mat'))
load(fullfile(DataPath(), 'ISIAmpDiscMultiFreqData_Processed.mat'))
load(fullfile(DataPath(), 'DetectionDataAll.mat'))

%% Reset
clearvars -except ISIAmpDiscData MechAmpDiscrimData DetectionData ISIAmpDiscMultiFreqData
%%% ISI 1 second data
% Make subject colors
num_ch = length(ISIAmpDiscData);
subj_colors = zeros(num_ch, 3);
for i = 1:num_ch
    subj_colors(i,:) = SubjectColors(ISIAmpDiscData(i).Subject);
end

[jnd1, pse1] = deal(zeros(num_ch, 3));
[djnd1, dpse1, pc2] = deal(zeros(num_ch,1));
for i = 1:num_ch
    jnd1(i,:) = ISIAmpDiscData(i).ISISigmoidSummary{"1","JND"};
    djnd1(i) = ISIAmpDiscData(i).ISISigmoidSummary{"1","dJND"};
    pse1(i,:) = ISIAmpDiscData(i).ISISigmoidSummary{"1","PSE"};
    dpse1(i) = ISIAmpDiscData(i).ISISigmoidSummary{"1","dPSE"} .* -1;
    pc2(i) = ISIAmpDiscData(i).ChoseI2;
end
dpse1_j = dpse1 ./ jnd1(:,1);

[mech_dpse1, mech_jnd] = deal(zeros(length(MechAmpDiscrimData), 1));
for i = 1:length(MechAmpDiscrimData)
    mech_dpse1(i) = MechAmpDiscrimData(i).SigmoidSummary{"1", "dPSE"} .* -1;
    mech_jnd(i) = MechAmpDiscrimData(i).SigmoidSummary{"1", "JND"}{1}(1);
end
mech_dpse1_j = mech_dpse1 ./ mech_jnd;

jnd_idx = jnd1(:, 1) < 40;
[p,h,s] = ranksum(jnd1(jnd_idx, 2), jnd1(jnd_idx, 3));
[p,h,s] = ranksum(pse1(jnd_idx, 2), pse1(jnd_idx, 3));
int_bias = abs(pc2 - 0.5);

[p,h,s] = ranksum(dpse1_j(jnd_idx), mech_dpse1_j);
[H,P,CI,STATS] = vartest2(dpse1_j(jnd_idx), mech_dpse1_j);

%%% Detection data
% Fix subject string
for i = 1:length(DetectionData)
    if strcmp(DetectionData(i).Subject, 'CRS02b')
        DetectionData(i).Subject = 'CRS02'; %#ok<SAGROW> 
    end
end

dts = NaN(num_ch, 1);
for i = 1:num_ch
    idx = find(strcmp({DetectionData.Subject}, ISIAmpDiscData(i).Subject) & [DetectionData.Channel] == ISIAmpDiscData(i).Channel);
    if ~isempty(idx)
        dts(i) = DetectionData(idx).MeanThreshold;
    end
end

%%% Frequency data
num_freq_ch = length(ISIAmpDiscMultiFreqData);
[freq_dpse, freq_djnd, freq_jnd] = deal(zeros(3, num_freq_ch));
for c = 1:num_freq_ch
    freq_dpse(:,c) = ISIAmpDiscMultiFreqData(c).Freq_PSEs.dPSE;
    freq_djnd(:,c) = ISIAmpDiscMultiFreqData(c).Freq_JNDs.dJND;
    freq_jnd(:,c) = ISIAmpDiscMultiFreqData(c).Freq_JNDs.JND(:,1);
end

[dpse_P,dpse_T,dpse_S] = anova2(freq_dpse);
[djnd_P,djnd_T,djnd_S] = anova2(freq_djnd);
[fjnd_P,fjnd_T,fjnd_S] = anova2(freq_jnd);

%% Main plot
% 2AFC task
% Psychometric function
% JND range
% Delta JND1
% Delta PSE
% Natural touth
SetFont('Arial', 9)

% Pulse config
px = [0, 100, 100, 300, 300, 400, 400, 800, 800, 900] ./ 1e6; % Phase durations
py = [0, 0, -1, -1, 0, 0, 0.5, 0.5, 0, 0]; % Relative amplitudes
np = 8; % Hz

clf; 
% 2AFC
axes('Position', [0.1 0.585 0.3 0.35]); hold on
    xl = [-1.25 3.25];
    yl = [-100 100];
    fc = 0.15;
    fcs = [range(xl)*fc, range(yl)*fc];

    % Trains
    px1 = repmat(px, [np,1])' + [0:np-1]./ np;
    py1 = repmat(py .* 60, [1, np]);
    py2 = repmat(py .* 80, [1, np]);
    plot([-1; px1(:); px1(:) + 2; 3], [0, py1, py2, 0], 'k')

    xlabel('Time (seconds)'); ylabel(sprintf('Amplitude (%sA)', GetUnicodeChar('mu')))

    % Crosses
    fixation_cross(-0.5, 100, fcs, 'w', 1.5)
    text(-0.5, 70, 'Pre', 'HorizontalAlignment','center','VerticalAlignment','middle')
    fixation_cross(0.5, 100, fcs, 'green', 1.5)
    text(0.5, 70, 'Int1', 'HorizontalAlignment','center','VerticalAlignment','middle')
    fixation_cross(1.5, 100, fcs, 'w', 1.5)
    text(1.5, 70, 'ISI', 'HorizontalAlignment','center','VerticalAlignment','middle')
    fixation_cross(2.5, 100, fcs, 'green', 1.5)
    text(2.5, 70, 'Int2', 'HorizontalAlignment','center','VerticalAlignment','middle')

    set(gca, 'YLim', yl, ...
             'XLim', xl, ...
             'XTick', [-1:3], ...
             'Clipping', 'off')

% Example psychometric function
e1 = 17;
axes('Position', [0.475 0.585 0.225 0.35]); hold on
    ordered_plot(ISIAmpDiscData(e1).ISISigmoidSummary, ISIAmpDiscData(e1).ConditionMeanTable, 1, gca())
    xlabel(sprintf('Comp. Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel('p(Comp. Higher)')
    set(gca, 'XTick', [40:10:80], ...
             'XLim', [35, 85])
    [x,y] = GetAxisPosition(gca, 5, 95);
    text(x, y, ColorText({'Std_1', 'All', 'Std_2'}, [.4 .4 1; .6 .6 .6; 1 .4 .4]), 'VerticalAlignment', 'top')

e2 = 12; %12
axes('Position', [0.725 0.585 0.225 0.35]); hold on
    ordered_plot(ISIAmpDiscData(e2).ISISigmoidSummary, ISIAmpDiscData(e2).ConditionMeanTable, 1, gca())
    xlabel(sprintf('Comp. Amplitude (%sA)', GetUnicodeChar('mu')))
    set(gca, 'XTick', [40:10:80], ...
             'XLim', [35, 85],...
             'YTickLabel', {})
    [x,y] = GetAxisPosition(gca, 5, 95);

% delta PSE range
axes('Position', [0.1 0.11 0.215 0.325]); hold on 
    plot([0 100], [0 100], 'Color', [.6 .6 .6], 'LineStyle', '--')
    pse_limit = [30, 90];
    pse_x = pse1(jnd_idx,2); pse_x(pse_x < pse_limit(1)) = pse_limit(1); pse_x(pse_x > pse_limit(2)) = pse_limit(2);
    pse_y = pse1(jnd_idx,3); pse_y(pse_y < pse_limit(1)) = pse_limit(1); pse_y(pse_y > pse_limit(2)) = pse_limit(2);
    scatter(pse_x, pse_y, 50, subj_colors(jnd_idx,:), 'filled')
    xlabel(sprintf('PSE_1 (%sA)', GetUnicodeChar('mu')), 'Color', [.4 .4 1])
    ylabel(sprintf('PSE_2 (%sA)', GetUnicodeChar('mu')), 'Color', [1 .4 .4])

    % add examples
    text(pse1(e1,2)+2, pse1(e1,3)+2, 'B', 'FontWeight', 'bold')
    text(pse1(e2,2)+2, pse1(e2,3)-2, 'C', 'FontWeight', 'bold')

    set(gca, 'XTick', [pse_limit(1):20:pse_limit(2)], ...
             'YTick', [pse_limit(1):20:pse_limit(2)], ...
             'XLim', pse_limit, ...
             'YLim', pse_limit)

% delta PSE/std ICMS vs Mech
axes('Position', [0.425 0.11 0.2 0.325]); hold on
    Swarm(1, dpse1_j, "Color", [.6 .6 .6], "SwarmColor", subj_colors, "DistributionStyle", "Box")
    Swarm(2, mech_dpse1_j, "Color", SubjectColors('BCI02'), "DistributionColor", [.6 .6 .6], "DistributionStyle", "Box",...
        "CenterColor", [.6 .6 .6])
    ylabel(sprintf('%sPSE / JND', GetUnicodeChar('Delta')))
    set(gca, 'XTick', [1,2], ...
             'XTickLabel', {'ICMS', 'Hand'}, ...
             'YTick', [-4:2:4], ...
             'XLim', [0.5 2.5], ...
             'YLim', [-4  4])
    text(1.5,-2, ColorText({'C1', 'P2', 'P3'}, [SubjectColors('BCI02'); SubjectColors('CRS02'); SubjectColors('CRS07')]), ...
        'HorizontalAlignment', 'center')

% |PSE| vs JND
axes('Position', [0.725 0.11 0.215 0.325]); hold on
    scatter(jnd1(jnd_idx,1), abs(dpse1(jnd_idx)), 50, subj_colors(jnd_idx,:), "filled")
    lm = fitlm(jnd1(jnd_idx,1), abs(dpse1(jnd_idx)), 'linear');
    x = linspace(5, 30)';
    y_pred = predict(lm, x);
    plot(x, y_pred, 'Color', [.4 .4 .4], 'LineStyle', '--', 'LineWidth', 2)
    xlabel(sprintf('JND (%sA)', GetUnicodeChar('mu')))
    ylabel(sprintf('|%sPSE| (%sA)', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))

set(gcf, 'Units', 'inches', ...
         'Position', [1, 1, 6.4, 4])

AddFigureLabels(gcf, [0.05, -0.015])
shg

%% Supplement plot
clf;
% JND range
axes('Position', [0.075 0.575 0.15 0.35]); hold on 
    Swarm(1, jnd1(:,1), 'SC', subj_colors, 'SwarmYLimits', [0 60], 'DS', 'Box')
    ylabel(sprintf('JND (%sA)', GetUnicodeChar('mu')))
    set(gca, 'XTick', [0], ...
             'XLim', [0.5 1.5])
    text(1.5, 60, ColorText({'C1', 'P2', 'P3'}, [SubjectColors('BCI02'); SubjectColors('CRS02'); SubjectColors('CRS07')]), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')

% dPSE vs PI2
axes('Position', [0.325 0.575 0.25 0.35]); hold on 
    scatter(pc2(jnd_idx), dpse1(jnd_idx), 50, subj_colors(jnd_idx,:), "filled")
    lm = fitlm(pc2(jnd_idx), dpse1(jnd_idx), 'linear');
    x = linspace(0.3, 0.7)';
    y_pred = predict(lm, x);
    plot(x, y_pred, 'Color', [.4 .4 .4], 'LineStyle', '--', 'LineWidth', 2)
    set(gca,'XLim', [.3 .7], ...
            'YTick', [-60:30:60], ...
            'XTick', [.3:.1:.7])
    xlabel('p(I_2)')
    ylabel(sprintf('%sPSE (%sA)', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))

% JND1 vs JND2
axes('Position', [0.675 0.575 0.25 0.35]); hold on
        plot([0 60], [0 60], 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(jnd1(jnd_idx,2), jnd1(jnd_idx,3), 50, subj_colors(jnd_idx,:), 'filled')
    xlabel(sprintf('JND_1 (%sA)', GetUnicodeChar('mu')), 'Color', [.4 .4 1])
    ylabel(sprintf('JND_2 (%sA)', GetUnicodeChar('mu')), 'Color', [1 .4 .4])

    set(gca, 'XTick', [0:20:60], ...
             'YTick', [0:20:60], ...
             'XLim', [0 60], ...
             'YLim', [0 60])
    [x,y] = GetAxisPosition(gca, 95, 50);

% dPSE Frequency
axes('Position', [0.075 0.1 0.25 0.35]); hold on
    for i = 1:size(freq_dpse, 2)
        plot([1:size(freq_dpse, 1)], freq_dpse(:,i), 'Color', [.6 .6 .6], 'LineStyle', '--')
        scatter([1:size(freq_dpse, 1)], freq_dpse(:,i), 50, SubjectColors('BCI02'), 'filled')
    end
    set(gca, 'XLim', [0.6 3.4], 'XTick', [1:3], 'XTickLabel', {'50', '100', '200'})
    xlabel('Frequency (Hz)')
    ylabel(sprintf('%sPSE  (%sA)', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))

% dJND Frequency
axes('Position', [0.385 0.1 0.25 0.35]); hold on
    for i = 1:size(freq_djnd, 2)
        plot([1:size(freq_djnd, 1)], freq_djnd(:,i), 'Color', [.6 .6 .6], 'LineStyle', '--')
        scatter([1:size(freq_djnd, 1)], freq_djnd(:,i), 50, SubjectColors('BCI02'), 'filled')
    end
    set(gca, 'XLim', [0.6 3.4], 'XTick', [1:3], 'XTickLabel', {'50', '100', '200'})
    xlabel('Frequency (Hz)')
    ylabel(sprintf('%sJND  (%sA)', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))

% JND Frequency
axes('Position', [0.7 0.1 0.25 0.35]); hold on
    for i = 1:size(freq_jnd, 2)
        plot([1:size(freq_jnd, 1)], freq_jnd(:,i), 'Color', [.6 .6 .6], 'LineStyle', '--')
        scatter([1:size(freq_jnd, 1)], freq_jnd(:,i), 50, SubjectColors('BCI02'), 'filled')
    end
    set(gca, 'XLim', [0.6 3.4], 'XTick', [1:3], 'XTickLabel', {'50', '100', '200'})
    xlabel('Frequency (Hz)')
    ylabel(sprintf('JND  (%sA)', GetUnicodeChar('mu')))

AddFigureLabels(gcf, [0.075, 0.01])
set(gcf, 'Units', 'inches', ...
         'Position', [1, 1, 6.4, 4])
shg


%% Helper functions
function ordered_plot(sig, cmt, isi_oi, ax)
    % Subsample 1s isi
    sig = sig{num2str(isi_oi), "Coeffs"};
    cmt = cmt(cmt.ISI == isi_oi, :);
    % Find standard
    std = mode(cmt.Amp1);
    % Get comparison
    ua = unique(cmt.Amp1);
    ua = ua(ua ~= std);
    % Separate by order
    i1_idx = cmt.Amp1 ~= std;
    y1 = cmt.CompHigher(i1_idx);
    y2 = cmt.CompHigher(~i1_idx);
    ym = mean([y1, y2], 2);

    xq = linspace(30, 90, 100);
    sf = GetSigmoid(2);
    scatter(ua, ym, 30, [.4 .4 .4], 'filled')
    plot(xq, sf(sig{1}(1,:), xq), 'Color', [.4 .4 .4], 'LineWidth', 2, 'Parent', ax)
    scatter(ua, y2, 30, [.4 .4 1], 'filled')
    plot(xq, sf(sig{1}(2,:), xq), 'Color', [.4 .4 1], 'LineWidth', 2, 'Parent', ax)
    scatter(ua, y1, 30, [1 .4 .4], 'filled')
    plot(xq, sf(sig{1}(3,:), xq), 'Color', [1 .4 .4], 'LineWidth', 2, 'Parent', ax)
end

function fixation_cross(x,y,s,c,lw)
    % x
    x1 = x - s(1)/2;
    x2 = x - s(1)/6;
    x3 = x + s(1)/6;
    x4 = x + s(1)/2;
    % y
    y1 = y - s(2)/2;
    y2 = y - s(2)/6;
    y3 = y + s(2)/6;
    y4 = y + s(2)/2;
    % Patch
    patch([x2, x3, x3, x4, x4, x3, x3, x2, x2, x1, x1, x2, x2], ...
          [y4, y4, y3, y3, y2, y2, y1, y1, y2, y2, y3, y3, y4], ...
          c, 'LineWidth', lw)
end
