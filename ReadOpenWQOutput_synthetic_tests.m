
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read project results: HDF5
outputFolder = '/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/synthTestT_results';

% Runs to process
model_all = {...'crhm'...
            ...'summa'...
            'mizuroute'...
            };
% test = 2; % 2_nrTrans_instS_PorMedia
% test = 4; % 4_nrTrans_contS_PorMedia
% test = 6; % 6_nrTrans_instS_PorMedia_linDecay
% test = 8; % 8_nrTrans_contS_PorMedia_linDecay
% test = 9, 10, 11, 11.1 12, 13 % chemistry
% test = 14; % test_debug
test = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extraction settings
extrData_flag = true;
debugMode_flag = true;

% Plot settings
if test >= 9
    plot_TimeX_ConcY_perElement_flag = true;
    plot_ConcX_ZProfileY_perTime_flag = false;
else % transport tests
    plot_TimeX_ConcY_perElement_flag = false;
    plot_ConcX_ZProfileY_perTime_flag = true;
end

openwq_noWaterConc_Val = -9999;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over models
for model_i = 1:numel(model_all)

    model_name = model_all{model_i};
    
    if model_name == "crhm"
        % crhm
        results_dir = '/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/code_crhm/case_studies/synthetic_tests';  
        comptName = 'SOIL_RECHR';
        
    elseif model_name == "summa"
        % summa
        results_dir = '/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/Summa-openWQ/case_studies/synthetic_tests';
        if (test >= 9); comptName = 'SCALARAQUIFER';
        else; comptName = 'ILAYERVOLFRACWAT_SOIL';
        end
     elseif model_name == "mizuroute"
        % summa
        results_dir = '/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/mizuRoute/case_studies/synthetic_tests';
        comptName = 'RIVER_NETWORK_REACHES';
        plot_TimeX_ConcY_perElement_flag = true;
        plot_ConcX_ZProfileY_perTime_flag = false;
    end
    
    test10_units = 'G';
    
    if model_name == "crhm" || model_name == "summa" 
        ix_request = 1;
        
        if test == 4 || test == 8 
            tStart_sec = seconds([...
                60 * 60 * 24 * 15,...
                60 * 60 * 24 * 70,...
                60 * 60 * 24 * 120,...
                ]);
        elseif test == 2 || test == 6 
            tStart_sec = seconds([...
                ...0,...
                60 * 60 * 24 * 55,...
                60 * 60 * 24 * 85,...
                60 * 60 * 24 * 120,...
                ]);
        end
        
    elseif model_name == "mizuroute" 
        ix_request = 6430;
        
        if test == 4 || test == 8 
        tStart_sec = seconds([...
                                60 * 60 * 24 * 100,...
                                60 * 60 * 24 * 400,...
                                60 * 60 * 24 * 800,...
                                ]);
        elseif test == 2 || test == 6 
            tStart_sec = seconds([...
                                ...0,...
                                60 * 60 * 24 * 55,...
                                60 * 60 * 24 * 85,...
                                60 * 60 * 24 * 120,...
                                ]);
        end
    end
    
    if plot_ConcX_ZProfileY_perTime_flag == true
        %tStart_sec = seconds([... % day
        %                    0,...
        %                    60 * 60 * 24,...
        %                    60 * 60 * 24 * 3,...
        %                    60 * 60 * 24 * 8,...
        %                    60 * 60 * 24 * 20,...
        %                    60 * 60 * 24 * 40,...
        %                    60 * 60 * 24 * 160,...
        %                    60 * 60 * 24 * 400,...
        %                    ]);
        if model_name == "crhm" || model_name == "summa" 
            requestProfileAPI = containers.Map(...
                            {'Profile', 'layer_m_interval', 'ReverseAxis_XY', 'TimeStamps'},...
                            {...
                                'z(x=1,y=1)',...
                                0.006,...
                                {true, true}, ...
                                datetime('1-Sep-2017 12:15:00') + tStart_sec...
                             });
        elseif model_name == "mizuroute" 
            requestProfileAPI = containers.Map(...
                            {'Profile', 'layer_m_interval', 'ReverseAxis_XY', 'TimeStamps'},...
                            {...
                                'x(y=1,z=1)',...
                                0.006,...
                                {true, true}, ...
                                datetime('1-Sep-2017 12:15:00') + tStart_sec...
                             });
        end
    end
    
    if test == 2

        Synthetic_test = '2_nrTrans_instS_PorMedia';
        
        if model_name == "summa" 

            extractElm_info = {...
                strcat(comptName,'@SPECIES_A#MG|L'),[repelem(1,100);
                                                     repelem(1,100);
                                                     1:100]';
                };
            
        elseif model_name == "mizuroute"
        
            extractElm_info = {...
                strcat(comptName,'@SPECIES_A#MG|L'),[[6567,6724,6893];...
                                                    repelem(1,3);...
                                                    repelem(1,3)]';};
        end

    elseif test == 4

        Synthetic_test = '4_nrTrans_contS_PorMedia';
        
        if model_name == "summa" 
            
            extractElm_info = {...
            strcat(comptName,'@SPECIES_A#MG|L'),[repelem(1,100);

                                                 repelem(1,100);
                                                 1:100]';
            };
                                             
        elseif model_name == "mizuroute"
        
            extractElm_info = {...
                strcat(comptName,'@SPECIES_A#MG|L'),[[6567,6724,6893];...
                                                    repelem(1,3);...
                                                    repelem(1,3)]';};
        end
        %extractElm_info = {...
        %    strcat(comptName,'@SPECIES_A#MG|L'),[1:100;
        %                                         repelem(1,100);
        %                                         repelem(1,100)]';
        %                                         };
    
    elseif test == 6

    Synthetic_test = '6_nrTrans_instS_PorMedia_linDecay';

    extractElm_info = {...
        strcat(comptName,'@SPECIES_A#MG|L'),[repelem(1,100);
                                             repelem(1,100);
                                             1:100]';
        };
        
    elseif test == 8

        Synthetic_test = '8_nrTrans_contS_PorMedia_linDecay';
        
        if model_name == "mizuroute" 
        
            extractElm_info = {...
                strcat(comptName,'@SPECIES_A#MG|L'),[[6567,6724,6893];...
                                                    repelem(1,3);...
                                                    repelem(1,3)]';
                strcat(comptName,'@time_track#MG|L'),[[6567,6724,6893];...
                                                    repelem(1,3);...
                                                    repelem(1,3)]';
                };
        else % summa and chrm
             extractElm_info = {...
        strcat(comptName,'@SPECIES_A#MG|L'),[repelem(1,100);
                                             repelem(1,100);
                                             1:100]';
        };
        end

    elseif test == 9
        Synthetic_test = '9_batch_singleSp_1storder';

        extractElm_info = {...
            strcat(comptName,'@SPECIES_A#MG'),[ix_request,1,1]
            };

    elseif test == 10
        Synthetic_test = '10_batch_singleSp_2ndorder';
        
        extractElm_info = {...
            strcat(comptName,'@SPECIES_A#', test10_units),[ix_request,1,1]
            };

    elseif test == 11
        Synthetic_test = '11_batch_2species';

        extractElm_info = {...
            strcat(comptName,'@SPECIES_A#MG'),[ix_request,1,1]
            strcat(comptName,'@SPECIES_B#MG'),[ix_request,1,1]
            };

    elseif test == 11.1
        Synthetic_test = '11_1_batch_3species';

        extractElm_info = {...
            strcat(comptName,'@SPECIES_A#MG'),[ix_request,1,1]
            strcat(comptName,'@SPECIES_B#MG'),[ix_request,1,1]
            strcat(comptName,'@SPECIES_C#MG'),[ix_request,1,1]
            };

    elseif test == 12
        Synthetic_test = '12_batch_nitrogencycle';

        extractElm_info = {...
            strcat(comptName,'@NREF#MG'),[ix_request,1,1]
            strcat(comptName,'@NLAB#MG'),[ix_request,1,1]
            strcat(comptName,'@DON#MG'),[ix_request,1,1]
            strcat(comptName,'@DIN#MG'),[ix_request,1,1]
            };    

    elseif test == 13
        Synthetic_test = '13_batch_oxygenBODcycle';

        extractElm_info = {...
            strcat(comptName,'@BOD#MG'),[ix_request,1,1]
            strcat(comptName,'@DEFICIT_OXYG#MG'),[ix_request,1,1]
            strcat(comptName,'@DO#MG'),[ix_request,1,1]
            };
        
    elseif test == 14
        Synthetic_test = 'test_debug';

        extractElm_info = {...
            strcat(comptName,'@SPECIES_A#MG'),[ix_request,1,1]
            strcat(comptName,'@TIME_TRACK#MG'),[ix_request,1,1]
            };
    end


    %% ================================================================================================
    % DON'T CHANGE BELOW ==============================================================================
    % ================================================================================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTARCT DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if extrData_flag == true

        if (model_name == "crhm" || model_name == "summa")
            openwq_readfuncs_dir = '../../build/source/openwq/openwq/supporting_scripts/Read_Outputs/';
        elseif (model_name == "mizuroute")
            openwq_readfuncs_dir = '../../mizuroute/route/build/openwq/openwq/supporting_scripts/Read_Outputs/';
        end
        addpath(openwq_readfuncs_dir)

        folderpath = strcat(results_dir, '/', Synthetic_test,'/', model_name,'/Output_OpenWQ/');

        output_openwq_tscollect_all = read_OpenWQ_outputs(...
            openwq_readfuncs_dir,...        % Fullpath for needed functions
            folderpath,...                  % Provide fullpath to directory where the HDF5 files are located
            extractElm_info,...             % If the flag above is 1, then provide details about the variables to plot
            'HDF5',...                      % Output format 
            debugMode_flag);                % Debug mode

        % Save Results
        save(strcat(outputFolder,'/',model_name, '/', Synthetic_test,'.mat'), 'output_openwq_tscollect_all');

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % TimeX_ConcY_perElement
    if plot_TimeX_ConcY_perElement_flag == true

        plot_OpenWQ_outputs_TimeX_ConcY_perElement(...
            output_openwq_tscollect_all,...
            openwq_noWaterConc_Val);  % No water concentration flag in openwq

    end

    % ConcX_ZProfileY_perTime
    if plot_ConcX_ZProfileY_perTime_flag == true

         plot_OpenWQ_outputs_ConcX_ZProfileY_perTime(...
            output_openwq_tscollect_all,...
            requestProfileAPI,...     % profile request, e.g., "z(x=1,y=1)"
            openwq_noWaterConc_Val);  % No water concentration flag in openwq

    end
        
end