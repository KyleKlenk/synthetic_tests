
clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change forcing file

% Runs to process

test_num = 1; % 1, 9, 10, 11, 11.1 12, 13

Forcing_paramList = {
    {'LWRadAtm', 350},...
    {'SWRadAtm', 0},...
    {'airpres', 101325},...
    {'airtemp', 283.16},...
    {'pptrate', 10},...
    {'spechum', 0.001},...
    {'windspd', 0},...
    };


newNC_sufix = '_new';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_name = 'summa';
DataType_2exam = 'forcing_data';

% Folder dir
if test_num == 1; folderTest = '1_conserv_instant_SW';
elseif test_num == 9; folderTest = '9_batch_singleSp_1storder';
elseif test_num == 10; folderTest = '10_batch_singleSp_2ndorder';
elseif test_num == 11; folderTest = '11_batch_2species';    
elseif test_num == 11.1; folderTest = '11_1_batch_3species';   
elseif test_num == 12; folderTest = '12_batch_nitrogencycle';   
elseif test_num == 13; folderTest = '13_batch_oxygenBODcycle';   
end

% Forcing file
if test_num == 1; nc_forcing = 'openWQ_syntheticTests_Transp_forcing.nc';
elseif test_num == 9; nc_forcing = 'openWQ_syntheticTests_BGC_forcing.nc';
elseif test_num == 10; nc_forcing = 'openWQ_syntheticTests_BGC_forcing.nc';
elseif test_num == 11; nc_forcing = 'openWQ_syntheticTests_BGC_forcing.nc';    
elseif test_num == 11.1; nc_forcing = 'openWQ_syntheticTests_BGC_forcing.nc';   
elseif test_num == 12; nc_forcing = 'openWQ_syntheticTests_BGC_forcing.nc';   
elseif test_num == 13; nc_forcing = 'openWQ_syntheticTests_BGC_forcing.nc';   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = strcat(folderTest,'/summa/summa/',DataType_2exam,'/');
nc_file = strcat(folder,nc_forcing);

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
newNCfile = strcat(folder,extractBefore(nc_forcing,'.nc'),...
    newNC_sufix,...
    '.nc');

copyfile(nc_file, newNCfile)

%ncdisp(newNCfile)

%%%%%%%%%%
% change variables
%%%%%%%%%%

% change time (same as crhm's simulation so that results can be compared)
% crhm_dateStart = [2017, 7, 28, 12, 15, 0];
% crhm_dateEnd = [2019, 12, 20, 13, 45, 0];
% summa_dateStart = crhm_dateStart;
% sartSec = etime(crhm_dateStart, summa_dateStart);
% sartDays = sartSec / (24 * 60 * 60);
% endSec = etime(crhm_dateEnd, summa_dateStart);
% endDays = endSec / (24 * 60 * 60);
% newTime = (sartDays: 1/ (4 * 24): endDays);
% ncwrite(newNCfile,'time',newTime);
% change data_step
%ncwrite(newNCfile,'data_step',900);

% change variables
for p = 1:numel(Forcing_paramList)
    
    varName = Forcing_paramList{p}{1};
    varNewNum = Forcing_paramList{p}{2};
    varNum = ncread(newNCfile,varName);
    
    newVar = repelem(varNewNum, numel(varNum));
    
    ncwrite(newNCfile,varName, newVar);
    
end


% Plot
figure   
numPanels_y = ceil(numel(Forcing_paramList)/2);
numPanels_x = ceil(numel(Forcing_paramList)/numPanels_y);

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

% Plot
for p = 1:numel(Forcing_paramList)

    varVals_old = ncread(nc_file,Forcing_paramList{p}{1});
    varVals_new = ncread(newNCfile,Forcing_paramList{p}{1});

    subplot(numPanels_x, numPanels_y, p)
    plot(time, varVals_old, 'linewidth', 2)
    hold on
    plot(time, varVals_new, 'linewidth', 2)
    xlabel('time')
    ylabel(Forcing_paramList{p}{1})
    legend('old','new')
    datetick('x','keeplimits','keepticks')
    grid on

end
