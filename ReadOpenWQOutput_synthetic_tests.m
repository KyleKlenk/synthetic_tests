
% Read project results: HDF5

% Tests possible = [9, 10, 11, 11.1 12, 13];
ExtrVisData = true;
test = 13;

outputFolder = "/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/synthTestT_results/";
model_all = {'crhm',...
            'summa'};
        
model_i = 2;

model_name = model_all{model_i};

if model_i == 1
    comptName = 'SOIL_RECHR';
elseif model_i ==2
    comptName = 'SCALARAQUIFER';
end

if test == 1
    
    Synthetic_test = '1_conserv_instant_SW';
    
    extractElm_info = {...
        'SWE@SPECIES_A#MG/L',[1,1,1];...
        ...'RUNOFF@SPECIES_A#MG/L',[1,1,1];...
        ...'SSR@SPECIES_A#MG/L',[1,1,1];...
        ...'SD@SPECIES_A#MG/L',[1,1,1];...
        'SOIL_RECHR@SPECIES_A#MG/L',[1,1,1];...
        %'SOIL_LOWER@SPECIES_A#MG/L',[1,1,1];...
        ...'SURFSOIL@SPECIES_A#MG/L',[1,1,1];...
        ...'GW@SPECIES_A#MG/L',[1,1,1];...
        'SWE@SPECIES_A#KG',[1,1,1];...
        %'RUNOFF@SPECIES_A#KG',[1,1,1];...
        %'SSR@SPECIES_A#KG',[1,1,1];...
        %'SD@SPECIES_A#KG',[1,1,1];...
        'SOIL_RECHR@SPECIES_A#KG',[1,1,1];...
        %'SOIL_LOWER@SPECIES_A#KG',[1,1,1];...
        %'SURFSOIL@SPECIES_A#KG',[1,1,1];...
        %'GW@SPECIES_A#KG',[1,1,1];...
        };

elseif test == 9
    Synthetic_test = '9_batch_singleSp_1storder';

    extractElm_info = {...
        strcat(comptName,'@SPECIES_A#MG'),[1,1,1]
        };
    
elseif test == 10
    Synthetic_test = '10_batch_singleSp_2ndorder';

    extractElm_info = {...
        strcat(comptName,'@SPECIES_A#G'),[1,1,1]
        };
    
elseif test == 11
    Synthetic_test = '11_batch_2species';

    extractElm_info = {...
        strcat(comptName,'@SPECIES_A#MG'),[1,1,1];...
        strcat(comptName,'@SPECIES_B#MG'),[1,1,1];...
        };
    
elseif test == 11.1
    Synthetic_test = '11_1_batch_3species';

    extractElm_info = {...
        strcat(comptName,'@SPECIES_A#MG'),[1,1,1];...
        strcat(comptName,'@SPECIES_B#MG'),[1,1,1];...
        strcat(comptName,'@SPECIES_C#MG'),[1,1,1];...
        };

elseif test == 12
    Synthetic_test = '12_batch_nitrogencycle';

    extractElm_info = {...
        strcat(comptName,'@Nref#MG'),[1,1,1];...
        strcat(comptName,'@Nlab#MG'),[1,1,1];...
        strcat(comptName,'@DON#MG'),[1,1,1];...
        strcat(comptName,'@DIN#MG'),[1,1,1]
        };    
    
elseif test == 13
    Synthetic_test = '13_batch_oxygenBODcycle';

    extractElm_info = {...
        strcat(comptName,'@BOD#MG'),[1,1,1];...
        strcat(comptName,'@DEFICIT_OXYG#MG'),[1,1,1];...
        strcat(comptName,'@DO#MG'),[1,1,1]
        };
end


%% ================================================================================================
    % DON'T CHANGE BELOW ==============================================================================
    % ================================================================================================
if ExtrVisData == true
   
    addpath("/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/Summa-openWQ/build/source/openwq/openwq/supporting_scripts/Read_Outputs")
    openwq_readfuncs_dir = "/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/6_Projects/1_GWF/2_WIP/code/Summa-openWQ/build/source/openwq/openwq/supporting_scripts/Read_Outputs";

    plot_elemt_flag = true;

    folderpath = strcat(Synthetic_test,'/', model_name,'/Output_OpenWQ/');

    output_openwq_tscollect_all = read_OpenWQ_outputs(...
        openwq_readfuncs_dir,...    % Fullpath for needed functions
        folderpath,...              % Provide fullpath to directory where the HDF5 files are located
        plot_elemt_flag,...         % Flag to specify if to plot variables
        extractElm_info,...         % If the flag above is 1, then provide details about the variables to plot
        'HDF5',...                  % Output format
        true);    % Debug mode

    % Save Results
    save(strcat(outputFolder,'/',model_name, '/', Synthetic_test,'.mat'), 'output_openwq_tscollect_all');
    
end