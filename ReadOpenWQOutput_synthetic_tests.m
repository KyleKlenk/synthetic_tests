
% Read project results: HDF5

% Tests possible = [9, 10, 11, 11.1 12, 13];
ExtrVisData = true;
test = 12;

outputFolder = "/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/code_crhm/apply/Case_Studies/synthetic_tests/99_analytical_solutions/openWQ_results";

if test == 9
    Synthetic_test = '9_batch_singleSp_1storder';

    extractElm_info = {...
        'SOIL_RECHR@SPECIES_A#KG',[1,1,1];...
        %'SOIL_RECHR@SPECIES_B#KG',[1,1,1; 2,1,1];...
        %'SOIL_RECHR@SPECIES_C#KG',[1,1,1; 2,1,1];...
        %'SOIL_RECHR@SPECIES_A#MG|L',[1,1,1; 2,1,1];...
        %'SOIL_RECHR@SPECIES_B#MG|L',[1,1,1; 2,1,1];...
        %'SOIL_RECHR@SPECIES_C#MG|L',[1,1,1; 2,1,1];...
        };
    
elseif test == 10
    Synthetic_test = '10_batch_singleSp_2ndorder';

    extractElm_info = {...
        'SOIL_RECHR@SPECIES_A#KG',[1,1,1]
        };
    
elseif test == 11
    Synthetic_test = '11_batch_2species';

    extractElm_info = {...
        'SOIL_RECHR@SPECIES_A#KG',[1,1,1];...
        'SOIL_RECHR@SPECIES_B#KG',[1,1,1];...
        };
    
elseif test == 11.1
    Synthetic_test = '11_1_batch_3species';

    extractElm_info = {...
        'SOIL_RECHR@SPECIES_A#KG',[1,1,1];...
        'SOIL_RECHR@SPECIES_B#KG',[1,1,1];...
        'SOIL_RECHR@SPECIES_C#KG',[1,1,1];...
        };

elseif test == 12
    Synthetic_test = '12_batch_nitrogencycle';

    extractElm_info = {...
        'SOIL_RECHR@Nref#KG',[1,1,1];...
        'SOIL_RECHR@Nlab#KG',[1,1,1];...
        'SOIL_RECHR@DON#KG',[1,1,1];...
        'SOIL_RECHR@DIN#KG',[1,1,1]
        };    
    
elseif test == 13
    Synthetic_test = '13_batch_oxygenBODcycle';

    extractElm_info = {...
        'SOIL_RECHR@BOD#KG',[1,1,1];...
        'SOIL_RECHR@DEFICIT_OXYG#KG',[1,1,1];...
        'SOIL_RECHR@DO#KG',[1,1,1]
        };
end


%% ================================================================================================
    % DON'T CHANGE BELOW ==============================================================================
    % ================================================================================================
if ExtrVisData == true
   
    addpath("/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/code_crhm/openwq/supporting_scripts/Read_Outputs")
    openwq_readfuncs_dir = "/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/code_crhm/openwq/supporting_scripts/Read_Outputs/";

    plot_elemt_flag = true;

    folderpath = strcat(Synthetic_test,'/Output_OpenWQ/');

    output_openwq_tscollect_all = read_OpenWQ_outputs(...
        openwq_readfuncs_dir,...    % Fullpath for needed functions
        folderpath,...              % Provide fullpath to directory where the HDF5 files are located
        plot_elemt_flag,...         % Flag to specify if to plot variables
        extractElm_info,...         % If the flag above is 1, then provide details about the variables to plot
        'HDF5',...                  % Output format
        true);    % Debug mode

    % Save Results
    save(strcat(outputFolder,'/',Synthetic_test,'.mat'), 'output_openwq_tscollect_all');
    
end