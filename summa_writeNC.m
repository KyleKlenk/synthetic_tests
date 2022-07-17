
clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Runs to process
model_name = 'summa';

tests_all = [1,...
    ...9, 10, 11, 11.1 12, 13...
    ];

DataType_2exam_all = {...
    ....'outputs',...
    'forcing'};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for 

test_i = 1;
test = tests_all(test_i);

DataType_2exam = DataType_2exam_all{d};
if strcmp(DataType_2exam, 'output')
    % Outputs
    folder = '9_batch_singleSp_1storder/summa/summa/output/';
    nc_files_all = {'openWQ_synthTests_timestep.nc',...
                'celia1990_G1-1_timestep.nc'};
elseif strcmp(DataType_2exam, 'forcing')        
    % Forcing
    folder = '9_batch_singleSp_1storder/summa/summa/forcing_data/';
    nc_files_all = {'celia1990_forcing.nc'};
end

nc_file_i = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc_file = strcat(folder,nc_files_all{nc_file_i});

% To get information about the nc file
ncinfo(nc_file)
% to display nc file
ncdisp(nc_file)
% to read a vriable 'var' exisiting in nc file
myvar_summa = ncread(nc_file,'data_step');
time_summa = ncread(nc_file,'time');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(time_summa, myvar_summa)
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change nc file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% duplicate original nc file
newNCfile = strcat(folder,'openWQ_syntheticTests_BGC_forcing.nc');
copyfile(nc_file, newNCfile)

% change attributes
%fileattrib(newNCfile)
%ncwriteatt(newNCfile,'hruId','Dimensions','lala');

ncdisp(newNCfile)

%%%%%%%%%%
% change variables
%%%%%%%%%%

% change time (same as crhm's simulation so that results can be compared)
crhm_dateStart = [2017, 7, 28, 12, 15, 0];
crhm_dateEnd = [2019, 12, 20, 13, 45, 0];
summa_dateStart = crhm_dateStart;
sartSec = etime(crhm_dateStart, summa_dateStart);
sartDays = sartSec / (24 * 60 * 60);
endSec = etime(crhm_dateEnd, summa_dateStart);
endDays = endSec / (24 * 60 * 60);
newTime = (sartDays: 1/ (4 * 24): endDays);
ncwrite(newNCfile,'time',newTime);

% change data_step
ncwrite(newNCfile,'data_step',900);

% change variables
ParamList = {
    'LWRadAtm',...
    'SWRadAtm',...
    'airpres',...
    'airtemp',...
    'pptrate',...
    'spechum',...
    'windspd'};

for p = 1:numel(ParamList)
    
    varNum = ncread(newNCfile,ParamList{p});
    newVar = repelem(varNum(1), numel(newTime));
    
    ncwrite(newNCfile,ParamList{p}, newVar);
    
end

ncwriteatt(newNCfile,'time','units', ['days since ',...
    datestr(summa_dateStart,'yyyy-mm-dd HH:MM:SS')... format: '1995-01-01 00:00:00'
    ])

% Plot
figure

numPanels_y = ceil(numel(ParamList)/2);
numPanels_x = ceil(numel(ParamList)/numPanels_y);
for p = 1:numel(ParamList)
    
    varVals = ncread(newNCfile,ParamList{p});
    
    subplot(numPanels_x, numPanels_y, p)
    plot(newTime, varVals)
    xlabel('time')
    ylabel(ParamList{p})
    title(ParamList{p})

end

ncdisp(newNCfile)