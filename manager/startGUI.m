function startGUI(targetDir, dataFormat)
% A function to call the ISETBioCSF GUI
%
% Syntax:
%   startGUI(targetDir, dataFormat)
%
% Description:
%    A function call to start the Manager's GUI
%
% Inputs:
%    targetDir  - String. The target directory for data.
%    dataFormat - Struct. A structure containing the required data
%                 formatting information.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    targetDir = uigetdir(targetDir);

    [varNamesLevel1, varValuesLevel1, varFormatsLevel1, ...
        varEditablesLevel1, varWidthsLevel1] = ...
        scanDirToRetrieveVarValues(targetDir, dataFormat);  

    hFig = figure(1);
    clf;
    set(hFig, 'Position', [10 10 1670 500], 'MenuBar', 'none', ...
        'Color', [1 1 1]);

    level1TableWidth = 0.98;
    handles.tableLevel1ExistingExps = uitable(...
        'parent', hFig, ...
        'Data', varValuesLevel1, ...
        'ColumnName', varNamesLevel1, ...
        'ColumnFormat', varFormatsLevel1, ...
        'ColumnEditable', varEditablesLevel1, ...
        'ColumnWidth', varWidthsLevel1, ...
        'FontSize', 12, ..., ...
        'FontName', 'Menlo', ...
        'Units', 'normalized', ...
        'Position', [0.01 0.40 level1TableWidth 0.5], ...
        'CellEditCallback', @(h, e) disp([e.Indices e.NewData]));

    varEditablesNewExpsLevel1 = repmat(true, [1 numel(varNamesLevel1)]);
    varEditablesNewExpsLevel1(end - 1) = false;
    varValuesLevel1NewExpsLevel1 = varValuesLevel1(1, :);
    varValuesLevel1NewExpsLevel1{end - 1} = false;
    varValuesLevel1NewExpsLevel1{end} = true;

    handles.tableLevel1NewExps = uitable(...
        'parent', hFig, ...
        'Data', varValuesLevel1NewExpsLevel1, ...
        'ColumnName', varNamesLevel1, ...
        'ColumnFormat', varFormatsLevel1, ...
        'ColumnEditable', varEditablesNewExpsLevel1, ...
        'ColumnWidth', varWidthsLevel1, ...
        'FontSize', 12, ...
        'FontName', 'Menlo', ...
        'Units', 'normalized', ...
        'Position', [0.01 0.02 level1TableWidth 0.28], ...
        'CellEditCallback', @(h, e) disp([e.Indices e.NewData]));

    handles.addExperimentButton = uicontrol(...
        'Style', 'Pushbutton', ... 
        'Units', 'normalized', ...
        'Position', [0.01 0.30 0.1 0.07], ... 
        'FontSize', 12, ...
        'FontName', 'Menlo', ...
        'String', 'Add new experiment', ...
        'Callback', {@addExperimentLevel1, handles.tableLevel1NewExps, ...
        handles.tableLevel1ExistingExps});

    handles.computeNewExperimentsButton = uicontrol(...
        'Style', 'Pushbutton', ... 
        'Units', 'normalized', ...
        'Position', [0.50 0.30 0.15 0.07], ... 
        'FontSize', 12, ...
        'FontName', 'Menlo', ...
        'String', 'Compute new experiments', ...
        'Callback', {@computeNewExperiments, ...
        handles.tableLevel1NewExps, handles.tableLevel1ExistingExps});

end

function varValues = scanTargetDirectory(h, e, table, targetDir)
% Scan the specified target directory
%
% Syntax:
%   varValues = scanTargetDirectory(h, e, table, targetDir)
%
% Description:
%    Scan the specified target directory for the provided table.
%
% Inputs:
%    h         - Handle. The relevant GUI handle.
%    e         - String. Descriptive string for generated exceptions.
%    table     - Table. The table to retrieve variable values from.
%    targetDir - String. The desired target directory.
%
% Outputs:
%    varValues - Struct. A structure containing the desired variable values
%
% Optional key/value pairs:
%    None.
%

end

function computeNewExperiments(h, e, newExpsTable, existingExpsTable)
% Compute new experiments
%
% Syntax:
%   computeNewExperiments(h, e, newExpsTable, existingExpsTable)
%
% Description:
%    Compute new experiments, pulling data from an existing table
%
% Inputs:
%    h                 - Handle. The relevant GUI handle.
%    e                 - String. Descriptive string of generated exceptions
%    newExpsTable      - Table. The table for storing new experiment data.
%    existingExpsTable - Table. The table containing the existing
%                        experiment data.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    fprintf('Compute new exps, and update new and existing tables\n');
end

function [varNames, varValues, varFormats, varEditables, varWidths] = ...
    scanDirToRetrieveVarValues(targetDir, dataFormat)
% Scan the specified target directory to retrieve variable values
%
% Syntax:
%   [varNames, varValues, varFormats, varEditables, varWidths] = ...
%       scanDirToRetrieveVarValues(targetDir, dataFormat)
%
% Description:
%    Scan the specified target directory to retrieve information about the
%    variables therein.
%
% Inputs:
%    targetDir    - String. The target directory (filepath).
%    dataFormat   - Struct. A structure containing relevant data
%                   formatting information.
%
% Outputs:
%    varNames     - Cell. A cell array containing the variable names.
%    varValues    - Cell. A cell array containing the variable values.
%    varFormats   - Cell. A cell array containing the types for each of the
%                   variables listed in varNames.
%    varEditables - Array. A boolean array indicating whether the variables
%                   are editable.
%    varWidths    - Cell. A cell array containing the variable's widths.
%
% Optional key/value pairs:
%    None.
%

    targetDirInfo = dir(targetDir);
    filesNum = numel(targetDirInfo);

    varNames = {};
    varFormats = {};
    for k = 1:numel(dataFormat.level1)
        d = dataFormat.level1{k};
        varNames{numel(varNames) + 1} = d{1};
        varFormats{numel(varFormats) + 1} = d{2};
    end
    originalVarNamesNum = numel(varNames);

    varWidths = {};
    for k = 1:numel(varNames)
        varWidths = cat(2, varWidths, 50);
    end

    varWidths{1} = 50;
    varWidths{2} = 50;
    varWidths{end-1} = 120;

    dataFormat.level2

    varNames2 = {};
    varFormats2 = {};
    for k = 1:numel(dataFormat.level2)
        d = dataFormat.level2{k};
        varNames{originalVarNamesNum + 3 + k} = d{1};
    end
    originalVarNames2Num = numel(dataFormat.level2)

    tableRow = 0;
    for k = 1:filesNum
        if (targetDirInfo(k).isdir) && ...
                (~strcmp(targetDirInfo(k).name, '.')) && ...
                (~strcmp(targetDirInfo(k).name, '..'))
            fprintf('name: %s, date: %s \n', ...
                targetDirInfo(k).name, targetDirInfo(k).date);
            tableRow = tableRow + 1;

            for varIndex = 1:originalVarNamesNum
                i1 = strfind(targetDirInfo(k).name, varNames{varIndex});
                if (varIndex < originalVarNamesNum)
                    i2 = strfind(targetDirInfo(k).name, ...
                        varNames{varIndex + 1});
                else
                    i2 = numel(targetDirInfo(k).name);
                end
                varValues{tableRow, varIndex} = ...
                    targetDirInfo(k).name(i1 + ...
                    numel(varNames{varIndex}):i2 - 1);
                if (tableRow == 1)
                     varEditables(varIndex) = false;
                end
            end
            % Add the file exists/needed fields
            if (tableRow == 1)
                varNames{originalVarNamesNum + 1} = 'responses exist';
                varNames{originalVarNamesNum + 2} = 'responses compute';
                varNames{originalVarNamesNum + 3} = 'visualize';
                varFormats{originalVarNamesNum + 1} = 'logical';
                varFormats{originalVarNamesNum + 2} = 'logical';
                varFormats{originalVarNamesNum + 3} = 'logical';
                varWidths{originalVarNamesNum + 1} = 110;
                varWidths{originalVarNamesNum + 2} = 110;
                varWidths{originalVarNamesNum + 3} = 110;
                varEditables(originalVarNamesNum + 1) = false;
                varEditables(originalVarNamesNum + 2) = true;
                varEditables(originalVarNamesNum + 3) = true;
            end
            varValues{tableRow, originalVarNamesNum + 1} = true;
            varValues{tableRow, originalVarNamesNum + 2} = false; 
            varValues{tableRow, originalVarNamesNum + 3} = false; 

            targetFileLevel2Index = 0;
            level2FilesNum = 2;

            for kk = 0:level2FilesNum - 1
                targetFileLevel2Index = targetFileLevel2Index + 1;

                for varIndex = 1:originalVarNames2Num
                    if (varIndex == 1)
                        if (kk == 0)
                            val = 'isomerizations';
                        else
                            val = 'photocurrents';
                        end
                    elseif (varIndex == 2)
                        val = 2;
                    elseif (varIndex == 3)
                        val = '5-fold';
                    elseif (varIndex == 4)
                        val = 60;
                    end

                    varValues{tableRow + kk, ...
                        originalVarNamesNum + 3 + varIndex} = val;

                    if (tableRow == 1)
                         varEditables(...
                             originalVarNamesNum + 3 + varIndex) = false;
                         varFormats{...
                             originalVarNamesNum + 3 + varIndex} = 'char';
                         varWidths{...
                             originalVarNamesNum + 3 + varIndex} = 90;
                         varNames{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 1} = ...
                             'classifier exists';
                         varNames{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 2} = ...
                             'classifier compute';
                         varNames{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 3} = 'visualize';
                         varFormats{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 1} = 'logical';
                         varFormats{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 2} = 'logical';
                         varFormats{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 3} = 'logical';
                         varWidths{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 1} = 110;
                         varWidths{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 2} = 110;
                         varWidths{originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 3} = 110;
                         varEditables(originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 1) = false;
                         varEditables(originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 2) = true;
                         varEditables(originalVarNamesNum + 3 + ...
                             originalVarNames2Num + 3) = true;
                    end

                    varValues{tableRow, originalVarNamesNum + 3 + ...
                        originalVarNames2Num + 1} = true;
                    varValues{tableRow, originalVarNamesNum + 3 + ...
                        originalVarNames2Num + 2} = false; 
                    varValues{tableRow, originalVarNamesNum + 3 + ...
                        originalVarNames2Num + 3} = false;
                end

            end
            tableRow = tableRow + level2FilesNum-1;

        end
    end

    for k = 1:numel(varNames)
        varNames{k} = sprintf(['<html><center /><font size= + 0>%s', ...
            '<br />&nbsp;</font></html>'], varNames{k});
    end

end

function addExperimentLevel1(h, e, table, existingExpsTable)
% Add a level 1 experiment
%
% Syntax:
%   addExperimentLevel1(h, e, table, existingExpsTable)
%
% Description:
%    Add a new top level experiment to an existing experiments table.
%
% Inputs:
%    h         - Handle. The relevant GUI handle.
%    e         - String. Descriptive string for generated exceptions.
%    table     - Table. The table to contain the new experiment
%    existingExpsTable
%              - Table. The existing experiments table.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    % get existing table data
    theData = get(table, 'Data'); 
    % replicate last row
    nRows = size(theData, 1);
    theLastData = theData(nRows, :);
    % Make the 'data file: | exists' field false
    theLastData{numel(theLastData) - 1} = false;
    theLastData{numel(theLastData)} = true;
    theData = cat(1, theData, theLastData);
    theEditables = repmat(true, [1 numel(theLastData)]);
    % update table
    set(table, 'Data', theData, 'ColumnEditable', theEditables);
end
