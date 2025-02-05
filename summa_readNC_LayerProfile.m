
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPACE ANALYSIS
% Output only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runs to process
model_name = 'summa'; % scriot only for summa (not applicable to CRHM)
summaNoDataFlag = -9999;

tests_all = [2,...
    ... 2, 4, 6, 8, 9, 10, 11, 11.1 12, 13...
    ];

timestep_to_print = [...
    60 * 60 * 24 * 15/900,...
    60 * 60 * 24 * 70/900,...
    60 * 60 * 24 * 120/900,...
    ];

DataType_2exam_all = {'output'};

% if choosing output, then select the parameters to plot
Output_paramList = {
    'pptrate',...
    'averageRoutedRunoff',...
    'scalarSWE',...
    'mLayerVolFracIce',...
    'mLayerVolFracLiq',...
    'mLayerMatricHead',...
    'scalarSnowfall',...
    'scalarSnowSublimation',...
    'scalarGroundEvaporation',...
    'mLayerTranspire',...
    'scalarThroughfallSnow',...
    'scalarThroughfallRain',...
    'scalarSnowDrainage',...
    'scalarInfiltration',...
    'iLayerLiqFluxSoil',...
    'mLayerBaseflow',...
    'scalarSoilBaseflow',...
    'scalarSoilDrainage',...    
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over tests
for test_i = 1:numel(tests_all)
    
    test_num = tests_all(test_i);
    
    % Folder dir
    if test_num == 2; folderTest = '2_nrTrans_instS_PorMedia';
    elseif test_num == 4; folderTest = '4_nrTrans_contS_PorMedia';
    %elseif test_num == 6; folderTest = '6_nrTrans_instS_PorMedia_linDecay';
    elseif test_num == 8; folderTest = '8_nrTrans_contS_PorMedia_linDecay';
    elseif test_num == 9; folderTest = '9_batch_singleSp_1storder';
    elseif test_num == 10; folderTest = '10_batch_singleSp_2ndorder';
    elseif test_num == 11; folderTest = '11_batch_2species';    
    elseif test_num == 11.1; folderTest = '11_1_batch_3species';   
    elseif test_num == 12; folderTest = '12_batch_nitrogencycle';   
    elseif test_num == 13; folderTest = '13_batch_oxygenBODcycle';   
    end
    
    % Output file
    if test_num == 2; nc_output = 'openWQ_synthTests_timestep.nc';
    elseif test_num == 4; nc_output = 'openWQ_synthTests_timestep.nc';
    %elseif test_num == 6; nc_output = 'openWQ_synthTests_timestep.nc';
    elseif test_num == 8; nc_output = 'openWQ_synthTests_timestep.nc';
    elseif test_num == 9; nc_output = 'openWQ_synthTests_timestep.nc';
    elseif test_num == 10; nc_output = 'openWQ_synthTests_timestep.nc';
    elseif test_num == 11; nc_output = 'openWQ_synthTests_timestep.nc';    
    elseif test_num == 11.1; nc_output = 'openWQ_synthTests_timestep.nc';   
    elseif test_num == 12; nc_output = 'openWQ_synthTests_timestep.nc';   
    elseif test_num == 13; nc_output = 'openWQ_synthTests_timestep.nc';   
    end
    
  
    for d = 1:numel(DataType_2exam_all)
        
        DataType_2exam = DataType_2exam_all{d};
        
        folder = strcat(folderTest,'/summa/summa/',DataType_2exam,'/');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        nc_file = strcat(folder,nc_output);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        paramList = Output_paramList;
        
        % Get time
        time_secs = ncread(nc_file,'time');
        timeAtr = ncreadatt(nc_file,'time','units');
        timeAtr_split = split(timeAtr,'since');
        timeAtr_units = strtrim(timeAtr_split{1});
        
        timeAtr_starDate = strtrim(timeAtr_split{2});
        timeAtr_starDate(strfind(timeAtr_starDate,' -'):end) = '';
        timeAtr_starDate = datevec(strtrim(timeAtr_starDate));
        
        if strcmp(timeAtr_units,'seconds')
            time = datetime(timeAtr_starDate) + seconds(time_secs);
        elseif strcmp(timeAtr_units,'days')
            time = datetime(timeAtr_starDate) + days(time_secs);
        end
        
         %Extract data 
        hbar = parfor_progressbar(...
        numel(paramList),...
        ['Extracting summa variables. File: ', nc_file]);
    
        Warning_text = cell(numel(paramList),1);
    
        for p = 1:numel(paramList)
            
            % update waitbar
            hbar.iterate(1);

            varVals_all = ncread(nc_file,paramList{p});
            size_varVals = size(varVals_all);
            
            if numel(size_varVals) > 2
                varVals = permute(varVals_all(1,:,:),[2 3 1]);
            else
                varVals = 'NOT_MULTIDIMENSIONAL';
                Warning_text{p} = strcat(...
                    '<WARNING> Variable not-multidimensional (space): ',...
                    paramList{p});
            end
            
            varVals_compile{p} = varVals;

        end
        close(hbar)
        
        % Find index of nont empty entries (so, the ones that have
        % multi-dimensional data
        EmptyIndex = find(~cellfun(@isempty,Warning_text));
        nonEmptyIndex = find(cellfun(@isempty,Warning_text));
        
        % Print to console the variables that are uni-dimensional
        Warning_textRelevant = Warning_text(EmptyIndex);
        for w = 1:numel(EmptyIndex)
            disp(Warning_textRelevant{w});
        end
         
        % Get only multi-dimensional data
        varVals_compile_MultDim = varVals_compile(nonEmptyIndex);
        paramList_multiDim = paramList(nonEmptyIndex);
        
        % Plot 
        % can't be paralellized
    
        figure('Name', nc_file)  
        
        numPanels_y = ceil(numel(varVals_compile_MultDim)/2);
        numPanels_x = ceil(numel(varVals_compile_MultDim)/numPanels_y);
        
        layer_seq = 1:1:numel(varVals_compile_MultDim{1}(:,1));
        
        hbar = parfor_progressbar(...
        numel(varVals_compile_MultDim),...
        ['Plotting summa variables. File: ', nc_file]);
        
        
        for p = 1:numel(varVals_compile_MultDim)
            
            % update waitbar
            hbar.iterate(1);
            
            % get var
            varVals = varVals_compile_MultDim{p};
            
            % Convert summa's noData flag to nan so it doesn't show up
            varVals(varVals == summaNoDataFlag) = NaN;
            
            % Sometimes it varies, so needs to be updated for every
            % variable
            layer_seq = 1:1:numel(varVals(:,1));
            legend_cell = cell(numel(timestep_to_print),1);
            
            for tstep = 1:numel(timestep_to_print)
                
                % Slice variable (at timestep_to_print)
                varVals_layer = varVals(:,timestep_to_print(tstep));

                % plot
                hold on
                subplot(numPanels_x, numPanels_y, p)
                plot(varVals_layer, layer_seq)
                set(gca, 'Ydir', 'reverse'); 
                xlabel('Var')
                ylabel('layer')
                grid on
                title(paramList_multiDim{p})
                legend_cell{tstep} = num2str(timestep_to_print(tstep));
                
            end
            
            legend(legend_cell)

        end
        close(hbar)

    end

end