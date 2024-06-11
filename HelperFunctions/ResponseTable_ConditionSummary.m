function summary_table = ResponseTable_ConditionSummary(input_table)
    % norm_table = ResponseTable_ConditionSummary(input_table)
    % Looks for unique combinations of all columns except the last and then
    % makes a summary, similar to ConditionMean but bigger
    input_array = table2array(input_table); % Convert to array for easy boolean operations
    unique_conditions = unique(input_array(:,1:end-1), 'rows'); % Unique combinations of parameters
    condition_mean = zeros(size(unique_conditions,1), 1);
    condition_values = cell(size(unique_conditions,1), 1);
    for c = 1:size(unique_conditions,1)
        % Find the rows of the input array where the parameters match the condition of interest
        condition_idx = all(input_array(:,1:end-1) == unique_conditions(c,:), 2);
        condition_mean(c) = mean(input_array(condition_idx,end), 'omitnan');
        condition_values{c} = input_array(condition_idx,end);
    end
    % Convert back to table and copy input table names
    summary_table = array2table(unique_conditions, 'VariableNames', input_table.Properties.VariableNames(1:end-1));
    summary_table.Mean = condition_mean;
    summary_table.Values = condition_values;
end
