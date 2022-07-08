
% Read CRHM outputs

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Directories
outputDir = "1_conserv_instant_SW/";
outputFile = "CRHM_output_1.txt";
obsFile = "weather_steady_state.obs";

% 2) Results
Var2Read_list = ["soil_rechr(1)", "soil_rechr(2)",...
                "soil_moist(1)","soil_moist(2)",...
                "SWE(1)","SWE(2)",...
                "snowmelt_int(1)", "snowmelt_int(2)"];
Var2Read_i = 7;

% 3) Observations (set to zero if unwanted)
Obs2Read_i = 6; % 6) Tair, 7) rh, 8) u, 9) p, 10) Qsi

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extraction code (Results)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Var2Read = Var2Read_list(Var2Read_i);

% Full path
outputFile_fullpath = strcat(outputDir, outputFile);

% Read file
fid = fopen(outputFile_fullpath);
hdr_res = strtrim(regexp(fgetl(fid),'\t','split'));
unt_res = strtrim(regexp(fgetl(fid),'\t','split'));
mat_res = cell2mat(textscan(fid,repmat('%f',1,numel(hdr_res))));
fclose(fid);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extraction code (Obs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Full path
obsFile_fullpath = strcat(outputDir, obsFile);

% Read file
fid = fopen(obsFile_fullpath);
frewind(fid);
fgetl(fid);
hdr_obs_1 = strtrim(regexp(fgetl(fid),'\t','split')); hdr_obs_1 = hdr_obs_1{1};
hdr_obs_2 = strtrim(regexp(fgetl(fid),'\t','split')); hdr_obs_2 = hdr_obs_2{1};
hdr_obs_3 = strtrim(regexp(fgetl(fid),'\t','split')); hdr_obs_3 = hdr_obs_3{1};
hdr_obs_4 = strtrim(regexp(fgetl(fid),'\t','split')); hdr_obs_4 = hdr_obs_4{1};
hdr_obs_5 = strtrim(regexp(fgetl(fid),'\t','split')); hdr_obs_5 = hdr_obs_5{1};
hdr_obs = {hdr_obs_1, hdr_obs_2, hdr_obs_3, hdr_obs_4, hdr_obs_5};
mat_obs = cell2mat(textscan(fid,repmat('%f',1, numel(hdr_obs) + 5),'headerlines',1));
fclose(fid);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Results
n_var = find(hdr_res == Var2Read);
time_res = datenum(mat_res(:,1)) + 693960;
var_res = mat_res(:,n_var);
TS_res = timeseries(var_res, time_res);

% Obs
var_obs = mat_obs(:,Obs2Read_i);
time_obs = datenum([mat_obs(:,1), mat_obs(:,2), mat_obs(:,3), mat_obs(:,4), mat_obs(:,5), zeros(numel(mat_obs(:,5)),1)]);
TS_obs = timeseries(var_obs, time_obs);
obs_name = hdr_obs{Obs2Read_i - 5};

figure
subplot(121)
plot(TS_res)
ylabel(Var2Read)
xlabel("Time")
title(Var2Read)
datetick('x','keepticks')
grid on

subplot(122)
plot(TS_obs)
ylabel(obs_name)
xlabel("Time")
title(obs_name)
datetick('x','keepticks')
grid on