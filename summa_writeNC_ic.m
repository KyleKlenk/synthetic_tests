
clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change initial conditions file

% Runs to process

test_num = 1; % 1, 9, 10, 11, 11.1 12, 13

IC_attributes = {
    {'scalarv', 1},...
    {'hru', 1},...
    {'ifcToto', 101},...
    {'midToto', 100},...
    {'midSoil', 100},...
};

IC_variables = {
    {'scalarCanopyTemp', 290},...
    {'dt_init', 10},...
    {'iLayerHeight', 0.1},... multidimensional
    {'mLayerVolFracLiq', 0.1},... multidimensional
    {'mLayerVolFracIce', 0},... multidimensional
    {'mLayerTemp', 285.1600},... multidimensional
    {'scalarSnowDepth', 0},...
    {'scalarCanopyLiq', 0},...
    {'scalarSfcMeltPond', 0},...
    {'nSoil', 100},...
    {'scalarAquiferStorage', 0},...
    {'mLayerDepth', 0.0060},... multidimensional
    {'scalarCanairTemp', 286},...
    {'nSnow', 0},...
    {'scalarSnowAlbedo', 0.8200},...
    {'mLayerMatricHead', -10},... multidimensional
    {'scalarSWE', 0},...
    {'scalarCanopyIce', 0},...
    };

%ncread(nc_file,"")

newNC_sufix = '_new';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_name = 'summa';

% Folder dir
if test_num == 1; folderTest = '1_conserv_instant_SW';
elseif test_num == 9; folderTest = '9_batch_singleSp_1storder';
elseif test_num == 10; folderTest = '10_batch_singleSp_2ndorder';
elseif test_num == 11; folderTest = '11_batch_2species';    
elseif test_num == 11.1; folderTest = '11_1_batch_3species';   
elseif test_num == 12; folderTest = '12_batch_nitrogencycle';   
elseif test_num == 13; folderTest = '13_batch_oxygenBODcycle';   
end

% Inicial conditions file
if test_num == 1; nc_initCond = 'summa_zInitialCond_OpenWQ_systheticTests_BGQ.nc';
elseif test_num == 9; nc_initCond = 'summa_zInitialCond_OpenWQ_systheticTests_BGQ.nc';
elseif test_num == 10; nc_initCond = 'summa_zInitialCond_OpenWQ_systheticTests_BGQ.nc';
elseif test_num == 11; nc_initCond = 'summa_zInitialCond_OpenWQ_systheticTests_BGQ.nc';    
elseif test_num == 11.1; nc_initCond = 'summa_zInitialCond_OpenWQ_systheticTests_BGQ.nc';   
elseif test_num == 12; nc_initCond = 'summa_zInitialCond_OpenWQ_systheticTests_BGQ.nc';   
elseif test_num == 13; nc_initCond = 'summa_zInitialCond_OpenWQ_systheticTests_BGQ.nc';   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = strcat(folderTest,'/summa/summa/SUMMA/');
nc_file = strcat(folder,nc_initCond);

% To get information about the nc file
% ncinfo(nc_file)
% to display nc file
% ncdisp(nc_file)

% to read a vriable 'var' exisiting in nc file
% myvar_summa = ncread(nc_file,'data_step');
% time_summa = ncread(nc_file,'time');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change nc file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% duplicate original nc file
newNCfile = strcat(folder,extractBefore(nc_initCond,'.nc'),...
    newNC_sufix,...
    '.nc');

copyfile(nc_file, newNCfile)

%ncdisp(newNCfile)
%ncread(nc_file,"scalarCanopyTemp")

% change attributes
%for p = 1:numel(IC_attributes)
    
%    attrName = IC_attributes{p}{1};
%    attrNewVal = IC_attributes{p}{2};
    
%    ncwriteatt(newNCfile, attrName, attrNewVal);
%    ncwrite()
%end

% change variables
for p = 1:numel(IC_variables)
    
    varName = IC_variables{p}{1};
    varNewNum = IC_variables{p}{2};
    varNum = ncread(newNCfile,varName);
    
    newVar = repelem(varNewNum, numel(varNum));
    
    ncwrite(newNCfile,varName, newVar);
    
end


% Plot
figure   
numPanels_y = ceil(numel(IC_variables)/2);
numPanels_x = ceil(numel(IC_variables)/numPanels_y);

% Plot
for p = 1:numel(IC_variables)

    varVals_old = ncread(nc_file,IC_variables{p}{1});
    varVals_new = ncread(newNCfile,IC_variables{p}{1});
    
    numLayer = numel(varVals_old);
    layerSeq = 1:1:numLayer;
    
    subplot(numPanels_x, numPanels_y, p)
    plot(varVals_old, layerSeq, 'linewidth', 2)
    hold on
    plot(varVals_new, layerSeq, 'linewidth', 2)
    ylabel('vertical layer')
    xlabel(IC_variables{p}{1})
    legend('old','new')
    %datetick('x','keeplimits','keepticks')
    grid on

end
