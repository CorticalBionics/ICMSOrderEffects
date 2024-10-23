%%% ISI Analysis
addpath(genpath("HelperFunctions"))
load(fullfile(DataPath(),  'ISIAmpDiscData_Processed.mat'))
load(fullfile(DataPath(),  'ISIMagEstData_Processed.mat'))
SetFont('Arial', 9)
subj_list = [ISIAmpDiscData.Subject];

%% Get summary 1s vs 5s dPSE
clearvars -except ISIAmpDiscData ISIMagEstData subj_list
%%% ISI 1 second data
% Make subject colors
num_ch = length(ISIAmpDiscData);
subj_colors = zeros(num_ch, 3);
for i = 1:num_ch
    subj_colors(i,:) = SubjectColors(ISIAmpDiscData(i).Subject);
end

[dpse1, dpse5, jnd1] = deal(zeros(num_ch, 1));
for i = 1:num_ch
    jnd1(i) = ISIAmpDiscData(i).ISISigmoidSummary{"1","JND"}(1);
    dpse1(i) = ISIAmpDiscData(i).ISISigmoidSummary{"1","dPSE"};
    dpse5(i) = ISIAmpDiscData(i).ISISigmoidSummary{"5","dPSE"};
end
jnd_idx = jnd1 > 40;
iidx = strcmp(subj_list, 'BCI02')';
dpse1(jnd_idx) = NaN;
dpse5(jnd_idx) = NaN;
[p,h,s] = signrank(abs(dpse5(iidx)), abs(dpse1(iidx)));
[a,b,c,d] = vartest2(dpse1(iidx), dpse5(iidx));

%% Other ISIs
idx = [5, 8, 17, 21];
isi = [0.5, 1, 2.5, 5];
dpse_all = NaN(length(idx), length(isi));
for i = 1:length(idx)
    ii = str2double(ISIAmpDiscData(idx(i)).ISISigmoidSummary.Properties.RowNames);
    for j = 1:length(ii)
        iidx = ii(j) == isi;
        dpse_all(i,iidx) = abs(ISIAmpDiscData(idx(i)).ISISigmoidSummary.dPSE(j));
    end
end

%% Mag est stuff
isis = [1, 5];
cond_amps = [10, 40, 80];
mag_effect = NaN(length(ISIMagEstData), 4, 2, 5); % Ch, cond_amp, ISI, test_amp
mag_tbl = cell(size(ISIMagEstData));
for i = 1:length(ISIMagEstData)
    % Anova stuff
    idx = [ISIMagEstData(i).ResponseTable.Ch1] == ISIMagEstData(i).TestChannel;
    mag_tbl{i} = ISIMagEstData(i).ResponseTable(idx, ["Amp1", "ISI", "Amp2", "Ch2", "NormResponse"]);
    % Plot stuff
    for j = 1:length(isis)
        idx = [ISIMagEstData(i).SummaryTable.Ch1] == 0 & ...
              [ISIMagEstData(i).SummaryTable.ISI] == isis(j);
        mag_effect(i,1,j,:) = ISIMagEstData(i).SummaryTable.Mean(idx);
        for c = 1:length(cond_amps)
            idx = [ISIMagEstData(i).SummaryTable.Ch1] == ISIMagEstData(i).TestChannel & ...
                  [ISIMagEstData(i).SummaryTable.ISI] == isis(j) & ...
                  [ISIMagEstData(i).SummaryTable.Amp1] == cond_amps(c);
            mag_effect(i,c+1,j,:) = ISIMagEstData(i).SummaryTable.Mean(idx);
        end
    end
end
% Divide by catch for effect size
catch_responses = squeeze(mean(mag_effect(:,1,:,:), 3));
mag_effect_norm = NaN(length(ISIMagEstData), 2, length(cond_amps));
for i = 1:length(ISIMagEstData)
    mag_effect_norm(i,1,:) = mean(squeeze(mag_effect(i, 2:end, 1, :)) ./ catch_responses(i,:), 2);
    mag_effect_norm(i,2,:) = mean(squeeze(mag_effect(i, 2:end, 2, :)) ./ catch_responses(i,:), 2);
end
% Run anova
mag_tbl = cat(1, mag_tbl{:});
[p,tbl,stats] = anovan(mag_tbl{:,"NormResponse"}, mag_tbl{:,["Amp1", "ISI", "Amp2", "Ch2"]},...
        'varnames', {'CondAmp', 'ISI', 'TestAmp', 'TestElec'}, 'display', 'on');

%% Plot
clearvars ax
clf;
set(gcf, 'Units', 'inches', ...
         'Position', [1, 1, 6.48, 3.5])
% dPSE 1v5
ax(1) = axes('Position', [0.075 0.5 0.225 0.4]); hold on
    plot([0 25], [0 25], 'Color', [.6 .6 .6], 'LineStyle', '--')
    scatter(abs(dpse1), abs(dpse5), 50, subj_colors, 'filled')
    xlabel(sprintf('|%sPSE| (%sA) @ 1s ISI', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))
    ylabel(sprintf('|%sPSE| (%sA) @ 5s ISI', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))
    set(gca, 'XLim', [0 50], 'YLim', [0 50], 'XTick', [0:10:50], 'YTick', [0:10:50])
    [x,y] = GetAxisPosition(gca, 95, 35);
    text(x,y, ColorText({'C1', 'P2', 'P3'}, [SubjectColors('BCI02'); SubjectColors('CRS02'); SubjectColors('CRS07')]), ...
        'HorizontalAlignment', 'right')

e = 3;
% 1s Psychometric
axes('Position', [0.1 0.75 0.085 0.15]); hold on
ordered_plot(ISIAmpDiscData(e).ISISigmoidSummary, ISIAmpDiscData(e).ConditionMeanTable, 1, gca())
    set(gca, 'XTick', [], 'YTick', [], ...
             'XLim', [35, 85])
    [x,y] = GetAxisPosition(gca, 50, -10);
    text(x, y, '1s ISI', 'VerticalAlignment', 'top', 'Color', [.6 .6 .6], 'HorizontalAlignment', 'center')
    
% 5s Psychometric
axes('Position', [0.2125 0.75 0.085 0.15]); hold on
    ordered_plot(ISIAmpDiscData(e).ISISigmoidSummary, ISIAmpDiscData(e).ConditionMeanTable, 5, gca())
    set(gca, 'XTick', [], 'YTick', [], ...
             'XLim', [35, 85])
    [x,y] = GetAxisPosition(gca, 50, -10);
    text(x, y, '5s ISI', 'VerticalAlignment', 'top', 'Color', [.6 .6 .6], 'HorizontalAlignment', 'center')

% Magnitude rating 
ax(2) = axes('Position', [0.415 0.5 0.225 0.4]); hold on
    e = 2;
    % Catch
    idx = [ISIMagEstData(e).SummaryTable.ISI] == 1 & ...
          [ISIMagEstData(e).SummaryTable.Ch1] == 0;
    c1 = ISIMagEstData(e).SummaryTable.Values(idx);
    idx = [ISIMagEstData(e).SummaryTable.ISI] == 5 & ...
          [ISIMagEstData(e).SummaryTable.Ch1] == 0;
    c2 = ISIMagEstData(e).SummaryTable.Values(idx);
    c3 = cellfun(@(c1,c2) [c1; c2], c1, c2, 'UniformOutput', false);
    % Combine catches
    AlphaLine([40:10:80], c3, [.6 .6 .6], 'ErrorType', 'SEM')
    % 1s
    idx = [ISIMagEstData(e).SummaryTable.ISI] == 1 & ...
          [ISIMagEstData(e).SummaryTable.Ch1] == ISIMagEstData(e).TestChannel & ...
          [ISIMagEstData(e).SummaryTable.Amp1] == 80;
    AlphaLine([40:10:80], ISIMagEstData(e).SummaryTable.Values(idx), SubjectColors('BCI02'), 'ErrorType', 'SEM')
    % 5s
    idx = [ISIMagEstData(e).SummaryTable.ISI] == 5 & ...
          [ISIMagEstData(e).SummaryTable.Ch1] == ISIMagEstData(e).TestChannel & ...
          [ISIMagEstData(e).SummaryTable.Amp1] == 80;
    AlphaLine([40:10:80], ISIMagEstData(e).SummaryTable.Values(idx), rgb(33, 150, 243), 'ErrorType', 'SEM')
    % Formatting
    xlabel(sprintf('Test Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel('Normalized Intensity')
    set(gca, 'XLim', [38, 82], 'YLim', [0.75 1.25], 'YTick', [0.8 : 0.1 : 1.2])
    [x,y] = GetAxisPosition(gca, 95, 5);
    text(x,y, ColorText({'1s ISI', '5s ISI', 'Catch'}, [SubjectColors('BCI02'); rgb(33, 150, 243); .6 .6 .6]), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')

% Magnitude conditioning magnitude plot
ax(3) = axes('Position', [0.75 0.5 0.225 0.4]); hold on
    plot([0.5 3.5], [1 1], 'Color', [.6 .6 .6], 'LineStyle', '--')
    for i = 1:length(cond_amps)
        Swarm(i, squeeze(mag_effect_norm(:,1,i)), SubjectColors('BCI02'), 'DS', 'Box', 'DM', 'None')
    end
    for i = 1:length(cond_amps)
        Swarm(i, squeeze(mag_effect_norm(:,2,i)), rgb(33, 150, 243), 'DS', 'Box', 'DM', 'None')
    end
    set(gca, 'XLim', [0.5 3.5], ...
             'YLim', [0.8 1.25], ...
             'XTickLabels', {'10', '40', '80'})
    xlabel(sprintf('Conditioning Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel('Relative Intensity')

% Many ISIs
e = 5; % 8, 17, 21
isis = ISIAmpDiscData(e).ISIs;
o = 0;
ax(4) = axes('Position', [(0.05 + o) 0.075 0.12 0.25]); hold on
for i = 1:length(isis)
    if i > 1
        axes('Position', [(0.05 + o) 0.075 0.12 0.25]); hold on
    end
    ordered_plot(ISIAmpDiscData(e).ISISigmoidSummary, ISIAmpDiscData(e).ConditionMeanTable, isis(i), gca())
    set(gca, 'XTick', [], 'YTick', [], ...
             'XLim', [35, 85])
    o = o + 0.14;
    [x,y] = GetAxisPosition(gca, 5, 100);
    text(x,y, sprintf('%0.1f s', isis(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Color', [.6 .6 .6])

    if i == 1
        ylabel('p(Comp Higher)')
    end
end

% dPSE at all ISIs
ax(5) = axes('Position', [0.825 0.075 0.15 0.25]); hold on
    for i = 1:size(dpse_all,1)
        plot([0.5, 1, 2.5, 5], dpse_all(i,:), 'Color', [.6 .6 .6], 'LineWidth', 1)
    end
    plot([0.5, 1, 2.5, 5], mean(dpse_all,1), 'Color', 'k', 'LineWidth', 2)
    x = xlabel('ISI (s)');
    x.Position(2) = x.Position(2) + 12;
    ylabel(sprintf('%sPSE (%sA)', GetUnicodeChar('Delta'), GetUnicodeChar('mu')))


AddFigureLabels(ax, [0.06, 0])
shg

%% Functions
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