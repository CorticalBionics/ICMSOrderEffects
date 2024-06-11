% ISIAmpDisc Analysis
addpath(genpath("HelperFunctions"))
load(fullfile(DataPath(), 'ISIMagEstData_Processed.mat'))
load(fullfile(DataPath(), 'ISIMagEstMultiDurData_Processed.mat'))
SetFont('Arial', 9)

%% Reset
clearvars -except ISIMagEstData ISIMagEstMultiDurData
num_ch = length(ISIMagEstData);

%% Filter to 1s ISI and repeat ANOVA
for i = 1:num_ch
    ISIMagEstData(i).ResponseTable = ISIMagEstData(i).ResponseTable([ISIMagEstData(i).ResponseTable.ISI] == 1, :); %#ok<*SAGROW> 
    % Create summary table
    anova_input = ISIMagEstData(i).ResponseTable(:, ["Ch1", "Amp1", "ISI", "Ch2", "Amp2", "NormResponse"]);
    ISIMagEstData(i).SummaryTable = ResponseTable_ConditionSummary(anova_input);
    % Run an ANOVA
    [p,tbl,stats] = anovan(anova_input{:,"NormResponse"}, anova_input{:,["Ch1", "Amp1", "Amp2"]},...
        'varnames', {'CondElec', 'CondAmp', 'TestAmp'}, 'display', 'off');
    ISIMagEstData(i).AnovaTable = tbl;
end

%% Effect of amplitude
amp_effect = NaN(num_ch, 4, 5); % Ch, cond_amp, test_amp
amp_cell = cell(size(ISIMagEstData));
for i = 1:num_ch
    % Pull out means for plots
    idx = [ISIMagEstData(i).SummaryTable.Ch1] == 0;
    amp_effect(i,1,:) = ISIMagEstData(i).SummaryTable.Mean(idx);
    ux = unique(ISIMagEstData(i).ResponseTable.Amp1([ISIMagEstData(i).ResponseTable.Ch1] == ISIMagEstData(i).TestChannel));
    for j = 1:length(ux)
        idx = [ISIMagEstData(i).SummaryTable.Ch1] == ISIMagEstData(i).TestChannel & ...
          [ISIMagEstData(i).SummaryTable.Amp1] == ux(j);
        amp_effect(i,j+1,:) = ISIMagEstData(i).SummaryTable.Mean(idx);
    end
    % Store table
    idx = [ISIMagEstData(i).ResponseTable.Ch1] == 0 | [ISIMagEstData(i).ResponseTable.Ch1] == ISIMagEstData(i).TestChannel;
    amp_cell{i} = ISIMagEstData(i).ResponseTable(idx, ["Amp1", "Ch2", "Amp2", "NormResponse"]);
end

% normalize by the 0th conditioning amp
amp_effect = amp_effect ./ amp_effect(:,1,:);
% average across test amplitudes
amp_effect = mean(amp_effect(:,2:end,:), 3);
% Anova
tbl = cat(1, amp_cell{:});
[p,tbl,stats] = anovan(tbl{:,"NormResponse"}, tbl{:,["Ch2", "Amp1", "Amp2"]},...
        'varnames', {'TestElec', 'CondAmp', 'TestAmp'}, 'display', 'on');

%% Effect of duration
valid_dur = [0, 0.1, 0.5, 1];
dur_effect = NaN(num_ch, length(valid_dur), 5); % Ch, cond_dur, test_amp
dur_cell = cell(length(valid_dur), 1);
for i = 1:num_ch
    % Pull out means for plots
    idx = [ISIMagEstMultiDurData(i).SummaryTable.Dur] == 0;
    dur_effect(i,1,:) = ISIMagEstMultiDurData(i).SummaryTable.Mean(idx);
    for j = 1:length(valid_dur(2:end))
        idx = [ISIMagEstMultiDurData(i).SummaryTable.Dur] == valid_dur(j+1);
        dur_effect(i,j+1,:) = ISIMagEstMultiDurData(i).SummaryTable.Mean(idx);
    end
    % Store table
    idx = ismember([ISIMagEstMultiDurData(i).ResponseTable.Dur], valid_dur);
    dur_cell{i} = ISIMagEstMultiDurData(i).ResponseTable(idx, ["Dur", "Ch", "Amp2", "NormResponse"]);
end

% normalize by the 0th conditioning amp
dur_effect = dur_effect ./ dur_effect(:,1,:);
% average across test amplitudes
dur_effect = mean(dur_effect(:,2:end,:), 3);
% Anova
tbl = cat(1, dur_cell{:});
[p,tbl,stats] = anovan(tbl{:,"NormResponse"}, tbl{:,["Dur", "Ch", "Amp2"]},...
        'varnames', {'Duration', 'TestElec', 'TestAmp'}, 'display', 'on');

%% Effect of stimulating location
loc_effect = NaN(num_ch, 3, 5); % Ch, cond_amp, test_amp
loc_cell = cell(size(ISIMagEstData));
for i = 1:num_ch
    % Catch
    idx = [ISIMagEstData(i).SummaryTable.Ch1] == 0;
    loc_effect(i,1,:) = ISIMagEstData(i).SummaryTable.Mean(idx);
    % Same
    idx = [ISIMagEstData(i).SummaryTable.Ch1] == ISIMagEstData(i).TestChannel & ...
          [ISIMagEstData(i).SummaryTable.Amp1] == 80;
    loc_effect(i,2,:) = ISIMagEstData(i).SummaryTable.Mean(idx);
    % Different
    idx = [ISIMagEstData(i).SummaryTable.Ch1] ~= ISIMagEstData(i).TestChannel & ...
          [ISIMagEstData(i).SummaryTable.Ch1] ~= 0 & ...
          [ISIMagEstData(i).SummaryTable.Amp1] == 80;
    loc_effect(i,3,:) = ISIMagEstData(i).SummaryTable.Mean(idx);
    % ANOVA tbl
    idx = ([ISIMagEstData(i).ResponseTable.Ch1] == ISIMagEstData(i).TestChannel & ...
           [ISIMagEstData(i).ResponseTable.Amp1] == 80) | ...
          ([ISIMagEstData(i).ResponseTable.Ch1] ~= ISIMagEstData(i).TestChannel & ...
           [ISIMagEstData(i).ResponseTable.Ch1] ~= 0 & ...
           [ISIMagEstData(i).ResponseTable.Amp1] == 80);
    loc_cell{i} = ISIMagEstData(i).ResponseTable(idx, ["Ch1", "Ch2", "Amp2", "NormResponse"]);
    loc_cell{i}.Same = loc_cell{i}.Ch1 == ISIMagEstData(i).TestChannel;
end

% normalize by the 0th conditioning amp
loc_effect = loc_effect ./ loc_effect(:,1,:);
% average across test amplitudes
loc_effect = mean(loc_effect(:,2:end,:), 3);
% Anova
tbl = cat(1, loc_cell{:});
[p,tbl,stats] = anovan(tbl{:,"NormResponse"}, tbl{:,["Same", "Amp2"]},...
        'varnames', {'CondElec', 'TestAmp'}, 'display', 'on');

%% Main plot
clf;
axes('Position', [0.2 0.81 0.7 0.175]); hold on
    i = 1;
    % Catch
    idx = [ISIMagEstData(i).SummaryTable.Ch1] == 0;
    xx = ISIMagEstData(i).SummaryTable.Amp2(idx);
    yy = ISIMagEstData(i).SummaryTable.Values(idx);
    AlphaLine(xx, yy, [.6 .6 .6], 'ErrorType', 'SEM')

    % 80 uA conditioning
    idx = [ISIMagEstData(i).SummaryTable.Ch1] == ISIMagEstData(i).TestChannel & ...
          [ISIMagEstData(i).SummaryTable.Amp1] == 80;
    xx = ISIMagEstData(i).SummaryTable.Amp2(idx);
    yy = ISIMagEstData(i).SummaryTable.Values(idx);
    AlphaLine(xx, yy, SubjectColors('BCI02'), 'ErrorType', 'SEM')

    xlabel(sprintf('Test Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel('Normalized Intensity')
    set(gca, 'XLim', [38, 82], 'YLim', [0.8 1.25], 'YTick', [0.8 : 0.1 : 1.2])
    [x,y] = GetAxisPosition(gca, 5, 95);
    text(x, y, ColorText({'Cond', 'Catch'}, [SubjectColors('BCI02'); .6 .6 .6]), ...
        'VerticalAlignment','top', 'HorizontalAlignment','left')

axes('Position', [0.2 0.56 0.7 0.175]); hold on
    plot([0.5 3.5], [1 1], 'Color', [.6 .6 .6], 'LineStyle','--')
    c = ColorGradient(rgb(225, 190, 231), rgb(123, 31, 162), size(amp_effect,2));
    for i = 1:size(amp_effect,2)
        Swarm(i, amp_effect(:,i), c(i,:), 'DM', 'None', 'DS', 'Box')
    end
    xlabel(sprintf('Conditioning Amplitude (%sA)', GetUnicodeChar('mu')))
    ylabel('Relative Intensity')
    set(gca, 'XTick', [1:3], ...
             'XLim', [0.5 3.5], ...
             'YLim', [0.8, 1.25], ...
             'XTickLabels', {'10', '40', '80'})

axes('Position', [0.2 0.31 0.7 0.175]); hold on
    plot([0.5 3.5], [1 1], 'Color', [.6 .6 .6], 'LineStyle','--')
    c = ColorGradient(rgb(225, 190, 231), rgb(123, 31, 162), size(dur_effect,2));
    for i = 1:size(dur_effect,2)
        Swarm(i, dur_effect(:,i), c(i,:), 'DM', 'None', 'DS', 'Box')
    end
    xlabel('Conditioning Duration (s)')
    ylabel('Relative Intensity')
    set(gca, 'XTick', [1:3], ...
             'XLim', [0.5 3.5], ...
             'YLim', [0.8, 1.25], ...
             'XTickLabels', valid_dur(2:end))

axes('Position', [0.2 0.06 0.7 0.175]); hold on
    plot([0.5 2.5], [1 1], 'Color', [.6 .6 .6], 'LineStyle','--')
    for i = 1:num_ch
        plot([1,2], loc_effect(i,:), 'Color', [.6 .6 .6], 'LineStyle', ':')
        scatter([1,2], loc_effect(i,:), 50, SubjectColors('BCI02'), 'filled')
    end
    xlabel('Conditioning Electrode Location')
    ylabel('Relative Intensity')
    set(gca, 'XTick', [1:2], ...
             'XLim', [0.5 2.5], ...
             'YLim', [0.8, 1.25], ...
             'XTickLabels', {'Same', 'Different'})

set(gcf, 'Units', 'inches', ...
         'Position', [1, 1, 2.25, 7])

AddFigureLabels(gcf, [0.15, 0.0275])
shg

