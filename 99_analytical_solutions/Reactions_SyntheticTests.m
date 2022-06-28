
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch Reactor Tests -> Reaction synthetic tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% Tests possible = [9, 10, 11, 11,1, 12, 13];
test = 10;

%%%%%%%%%%%%%%%%%%%%%%%
% Test 9
%%%%%%%%%%%%%%%%%%%%%%%
if (test == 9)

    c0 = 10;
    k = 0.01;
    
    % OpenWQ
    openWQres_file = "openWQ_results/9_batch_singleSp_1storder.mat";
    conc_A_openwq = getOpenWQResults_1species(openWQres_file);
    
    % Analytical
    tsim = numel(conc_A_openwq.Time);
    conc_A_analytical = test_1_singleSpecies_1order(c0, k, tsim);
    
    title_str = "Batch Reactor: Single species with 1^{st} order decay";
    plotAnalytical_singleSpec(conc_A_analytical,...
                              conc_A_openwq,...
                              title_str)

%%%%%%%%%%%%%%%%%%%%%%%
% Test 10
%%%%%%%%%%%%%%%%%%%%%%%
elseif (test == 10)
    
    c0 = 10;
    k = 0.01;

    % Numerical
    openWQres_file = "openWQ_results/10_batch_singleSp_2ndorder.mat";
    conc_A_openwq = getOpenWQResults_1species(openWQres_file);
    
    % Analytical
    tsim = numel(conc_A_openwq.Time);
    conc_A_analytical = test_2_singleSpecies_2ndorder(c0, k, tsim);
    
    title_str = "Batch Reactor: Single species with 2^{nd} order decay";
    plotAnalytical_singleSpec(conc_A_analytical,...
                              conc_A_openwq,...
                              title_str)

%%%%%%%%%%%%%%%%%%%%%%%
% Test 11
%%%%%%%%%%%%%%%%%%%%%%%
elseif (test == 11)
    
    c0_a = 10;
    c0_b = 0;
    k_a = 0.03;
    k_b = 0.01;
    
    openWQres_file = "openWQ_results/11_batch_2species.mat";
    [conc_A_openwq, conc_B_openwq] = getOpenWQResults_2species(openWQres_file);
    
    % Analytical
    tsim = numel(conc_A_openwq.Time);
    [conc_A_analytical, conc_B_analytical] = test_3_twoSpecies(c0_a, c0_b, k_a, k_b, tsim);
    
    title_str = "Batch Reactor: Chain reaction with 2 Species with 1^{st} order decay";
    plotAnalytical_2species(conc_A_analytical,...
                            conc_B_analytical,...
                            conc_A_openwq,...
                            conc_B_openwq,...
                            title_str,...
                            "Species\_A",...
                            "Species\_B")

%%%%%%%%%%%%%%%%%%%%%%%
% Test 11.1
%%%%%%%%%%%%%%%%%%%%%%%
elseif (test == 11.1)
    
    c0_a = 10;
    c0_b = 0;
    c0_c = 0;
    k_a = 0.03;
    k_b = 0.01;
    k_c = 0.005;
    
    openWQres_file = "openWQ_results/11_1_batch_3species.mat";
    [conc_A_openwq, conc_B_openwq, conc_C_openwq] = getOpenWQResults_3species(openWQres_file);
    
    % Analytical
    tsim = numel(conc_A_openwq.Time);
    [conc_A_analytical, conc_B_analytical, conc_C_analytical] = test_3_threeSpecies(c0_a, c0_b, c0_c, k_a, k_b, k_c, tsim);
    
    title_str = "Batch Reactor: Chain reaction with 3 Species with 1^{st} order decay";
    plotAnalytical_3species(conc_A_analytical,...
                            conc_B_analytical,...
                            conc_C_analytical,...
                            conc_A_openwq,...
                            conc_B_openwq,...
                            conc_C_openwq,...
                            title_str,...
                            "Species\_A",...
                            "Species\_B",...
                            "Species\_C")
    
    
%%%%%%%%%%%%%%%%%%%%%%%
% Test 12
%%%%%%%%%%%%%%%%%%%%%%%
elseif (test == 12)

    c0_Nref = 10;
    c0_Nlab = 10;
    c0_DON = 2;
    c0_DIN = 5;
    c0_N2 = 0;
    c0_plants = 0;
    k_degrad = 0.006;
    dissol_1 = 0.0002;
    dissol_2 = 0.0003;
    miner = 0.003;
    denit = 0.001;
    plantup = 0.001;
    
    openWQres_file = "openWQ_results/12_batch_nitrogencycle.mat";
    [conc_Nref_openwq, conc_Nlab_openwq, conc_DON_openwq, conc_DIN_openwq] = getOpenWQResults_Ncyclespecies(openWQres_file);
    
    % Analytical
    tsim = numel(conc_Nref_openwq.Time);
    [conc_Nref_analytical, conc_Nlab_analytical, conc_DON_analytical, conc_DIN_analytical] = test_4_nitrogen(c0_Nref,...
        c0_Nlab,...
        c0_DON,...
        c0_DIN,...
        c0_N2,...
        c0_plants,...
        k_degrad,...
        dissol_1,...
        dissol_2,...
        miner,...
        denit,...
        plantup,...
        tsim);

    title_str = "Batch Reactor: Nitrogen Cycle";
    plotAnalytical_NcycleRes(conc_Nref_analytical,...
                            conc_Nlab_analytical,...
                            conc_DON_analytical,...
                            conc_DIN_analytical,...
                            conc_Nref_openwq,...
                            conc_Nlab_openwq,...
                            conc_DON_openwq,...
                            conc_DIN_openwq,...
                            title_str,...
                            "Nref",...
                            "Nlab",...
                            "DON",...
                            "DIN")
    
%%%%%%%%%%%%%%%%%%%%%%%
% Test 13
%%%%%%%%%%%%%%%%%%%%%%%
elseif (test == 13)
    
    c0_bod = 10;
    c0_do_deficit = 0;
    c0_do_sat = 12;
    k_bod = 0.1;
    k_rear = 0.5;
    
    % Model works with mass balance, we either work with state-variables (DO) or
    % anti-state-variables (BOD, DO_defitic). In this case we are using variables
    
    %%%%%%%%%%%%%%%%%%%%%
    % Anti-State-Variables
    %%%%%%%%%%%%%%%%%%%%
    
    % Numerical
    openWQres_file = "openWQ_results/13_batch_oxygenBODcycle.mat";
    [c_bod_openwq, c_do_openwq_deficit] = getOpenWQResults_2species(openWQres_file);
    
    % Analytical
    tsim = numel(c_bod_openwq.Time);
    [c_bod_analytical, c_do_deficit_analytical] = test_Streeter_Phelps(c0_bod, c0_do_deficit, k_bod, k_rear, tsim);
    
    
    %%%%%%%%%%%%%%%%%%%%%
    % Calculatinf State-Variables
    %%%%%%%%%%%%%%%%%%%%
    
    % Post calculating the relevant state-variable DO (opossite to the
    % anti-state variable used DO_deficit)
    c_do_analytical = c0_do_sat - c_do_deficit_analytical;
    c_do_openwq = c_do_openwq_deficit;
    c_do_openwq.data_save_final = c0_do_sat - c_do_openwq_deficit.data_save_final;
    
    title_str = "Batch Reactor: Streeter Phelps (DO sag)";
    plotAnalytical_2species(c_bod_analytical,...
                            c_do_analytical,...
                            c_bod_openwq,...
                            c_do_openwq,...
                            title_str,...
                            "BOD",...
                            "DO")
                        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = getOpenWQResults_1species(openWQres_file)

    openwq_data = load(openWQres_file);
    
    c = openwq_data.output_openwq_tscollect_all{1,2}{1,2};

end

function [ca,cb] = getOpenWQResults_2species(openWQres_file)

    openwq_data = load(openWQres_file);
    
    ca = openwq_data.output_openwq_tscollect_all{1,2}{1,2};
    cb = openwq_data.output_openwq_tscollect_all{1,2}{2,2};

end

function [ca,cb, cc] = getOpenWQResults_3species(openWQres_file)

    openwq_data = load(openWQres_file);
    
    ca = openwq_data.output_openwq_tscollect_all{1,2}{1,2};
    cb = openwq_data.output_openwq_tscollect_all{1,2}{2,2};
    cc = openwq_data.output_openwq_tscollect_all{1,2}{3,2};

end

function [conc_Nref_openwq, conc_Nlab_openwq, conc_DON_openwq, conc_DIN_openwq] = getOpenWQResults_Ncyclespecies(openWQres_file)
    
    openwq_data = load(openWQres_file);

    conc_Nref_openwq = openwq_data.output_openwq_tscollect_all{1,2}{1,2};
    conc_Nlab_openwq = openwq_data.output_openwq_tscollect_all{1,2}{2,2};
    conc_DON_openwq = openwq_data.output_openwq_tscollect_all{1,2}{3,2};
    conc_DIN_openwq = openwq_data.output_openwq_tscollect_all{1,2}{4,2};

end

% Single species with 1st order decay
function c = test_1_singleSpecies_1order(c0, k, tsim)
    
    c = zeros(tsim, 1);
    
    for t = 1:tsim
        ti = t - 1;
        c(t) = c0 * exp(-k * ti);
    end
    
    
end

% Single species with second order decay
function c = test_2_singleSpecies_2ndorder(c0, k, tsim)
    
    c = zeros(tsim, 1);
    
    for t = 1:tsim
        ti = t - 1;
        c(t) = 1 / (1/c0 + k * ti);
    end
    
    
end


function [c_a, c_b] = test_3_twoSpecies(c0_a, c0_b, k_a, k_b, tsim)

    c_a = zeros(tsim, 1);
    c_b = zeros(tsim, 1);
    
    c_a_numerical = zeros(tsim, 1);
    c_a_numerical(1) = c0_a;
    c_b_numerical = zeros(tsim, 1);
    c_b_numerical(1) = c0_b;
    
    % Species 1 (simple single species linear decay)
    for t = 1:tsim
        
        ti = t - 1;
        
        % Species 1
        c_a(t) = c0_a * exp(-k_a * ti);
        
        % Species 2
        A = k_a * c0_a / (k_b - k_a);
        B = exp(-k_a * ti) - exp(-k_b * ti);
        C = c0_b * exp(-k_b * ti);
        
        c_b(t) = A * B + C;
        
        % Numerical approximation for testing purposes
        dca_dt = k_a * c_a_numerical(t) * (t-ti);
        c_a_numerical(t+1) = c_a_numerical(t) - dca_dt;
        dcb_dt = k_b * c_b_numerical(t) * (t-ti);
        c_b_numerical(t+1) = c_b_numerical(t) + dca_dt - dcb_dt;
        
    end
    
end

function [c_a, c_b, c_c] = test_3_threeSpecies(c0_a, c0_b, c0_c, k_a, k_b, k_c, tsim)
    
    c_a = zeros(tsim, 1);
    c_b = zeros(tsim, 1);
    c_c = zeros(tsim, 1);
    
    % Species 1 (simple single species linear decay)
    for t = 1:tsim
        
        ti = t - 1;
        
        % Species 1
        c_a(t) = c0_a * exp(-k_a * ti);
        
        % Species 2
        A = k_a * c0_a / (k_b - k_a);
        B = exp(-k_a * ti) - exp(-k_b * ti);
        C = c0_b * exp(-k_b * ti);
        
        c_b(t) = A * B + C;
        
        %Species 3
        A = (k_a * c0_a) / (k_b  - k_a);
        B = exp(-k_a * t) / (k_c - k_a);
        C = exp(-k_b * t) / (k_c - k_b);
        D = (c0_b * exp(-k_b * t)) / (k_c - k_b);
        E = c0_c * exp(-k_c * t);
        F = exp(-k_c * t) / (k_c - k_a);
        G = exp(-k_c * t) / (k_c - k_b);
        
        c_c(t) = k_b * (A * (B - C) + D) + E - k_b * (A * (F - G) + D);
        
    end

    %figure
    %plot(c_a, 'k-', 'linewidth',1.5)
    %hold on
    %plot(c_b, 'k--', 'linewidth',1.5)
    %hold on
    %plot(c_c, 'k:', 'linewidth',1.5)
    %xlabel("Time")
    %ylabel("Concentration")
    %legend("Species A",...
    %        "Species B",...
    %        "Species C")
    %grid on
    
    
end

% Nitrogen Cycle
function [conc_Nref_analytical, conc_Nlab_analytical, conc_DON_analytical, conc_DIN_analytical] = ...
    test_4_nitrogen(...
    c0_Nref,...
    c0_Nlab,...
    c0_DON,...
    c0_DIN,...
    c0_N2,...
    c0_plants,...
    k_degrad,...
    dissol_1,...
    dissol_2,...
    miner,...
    denit,...
    plantup,...
    tsim)
    
    
    conc_Nref_analytical = zeros(tsim, 1);
    conc_Nlab_analytical = zeros(tsim, 1);
    conc_DON_analytical = zeros(tsim, 1);
    conc_DIN_analytical = zeros(tsim, 1);
    
    conc_Nref_numerical = zeros(tsim, 1);
    conc_Nlab_numerical = zeros(tsim, 1);
    conc_DIN_numerical = zeros(tsim, 1);
    conc_Nref_numerical(1) = c0_Nref;
    conc_Nlab_numerical(1) = c0_Nlab;
    conc_DIN_numerical(1) = c0_DIN;
    
    % Species 1 (simple single species linear decay)
    for t = 1:tsim
        
        ti = t - 1;
        
        % Nref
        c0_a = c0_Nref;
        k_a = k_degrad + dissol_1;
  
        conc_Nref_analytical(t) = c0_a * exp(-k_a * ti);
        
        % Nlab
        c0_a = c0_Nref;
        c0_b = c0_Nlab;
        k_a = k_degrad;
        k_b = dissol_2 + miner;
        
        A = k_a * c0_a / (k_b - k_a);
        B = exp(-k_a * ti) - exp(-k_b * ti);
        C = c0_b * exp(-k_b * ti);
        conc_Nlab_analytical(t) = A * B + C;

        c0_a = c0_Nref;
        c0_b = c0_Nlab;
        c0_c = c0_DON;
        k_a1 = k_degrad;
        k_a2 = dissol_1;
        k_b1 = miner;
        k_b2 = dissol_2;
        
        A = ( (k_a2 * c0_a) / (-( k_a1 + k_a2 )) ) * exp( - (k_a1 + k_a2) * ti);
        B = (k_b2 * k_a1 * c0_a) / (k_b1 + k_b2 - k_a1) * (exp(-k_a1 * ti) / (-k_a1) - exp( -(k_b1 + k_b2) * ti) / (-(k_b1 + k_b2)));
        C = ( k_b2 * c0_b ) / (-(k_b1 + k_b2)) * exp(-(k_b1 + k_b2) * ti);
        
        const1 = (k_a2 * c0_a) / (k_a1 + k_a2);
        const2 = (k_b2 * k_a1 * c0_a) / (k_b1 + k_b2 - k_a1) * ( (1/(k_b1 + k_b2)) - (1 / k_a1) );
        const3 = k_b2 * c0_b / (k_b1 + k_b2);
        const = c0_c + const1 - const2 + const3;
        
        conc_DON_analytical(t) = A + B + C + const;
        
        %DIN (as the 3 species of the 3-species chain)
        %c0_a = c0_Nref;
        %c0_b = c0_Nlab;
        %c0_c = c0_DIN;
        %k_a = k_degrad;
        %k_b = miner + dissol_2;
        %k_c = denit + plantup;
        
        %A = (k_a * c0_a) / (k_b  - k_a);
        %B = exp(-k_a * t) / (k_c - k_a);
        %C = exp(-k_b * t) / (k_c - k_b);
        %D = (c0_b * exp(-k_b * t)) / (k_c - k_b);
        %E = c0_c * exp(-k_c * t);
        %F = exp(-k_c * t) / (k_c - k_a);
        %G = exp(-k_c * t) / (k_c - k_b);
        %conc_DIN_analytical(t) = k_b * (A * (B - C) + D) + E - k_b * (A * (F - G) + D);
        
        % New derivation
        cx0 = c0_Nref;
        A0 = c0_Nlab;
        ka = miner;
        kb = denit + plantup;
        k1 = k_degrad;
        k2 = miner;
        k3 = dissol_2;
        
        A = (ka * k1 * cx0) / (k2 + k3 - k1);
        B = (exp(-k1 * ti)) / (kb - k1) - (exp((-k2-k3)*ti)) / (kb - k2 - k3);
        C = (ka * A0) / (kb - k2 - k3);
        D = exp((-k2 - k3) * ti);
        
        const1 = A;
        const2 = 1 / (kb - k1) - 1 / (kb - k2 - k3);
        const3 = C;
        const = c0_DIN - const1 * const2 - const3;
        E = const * exp(-kb * t);
        
        conc_DIN_analytical(t) = A * B + C * D + E;
        
        % DIN numerical
        %if t < tsim
        %    conc_Nref_numerical(t + 1) = conc_Nref_numerical(t)...
        %        - conc_Nref_numerical(t) * k_degrad;
        %    conc_Nlab_numerical(t + 1) = conc_Nlab_numerical(t) ...
        %        + conc_Nref_numerical(t) * k_degrad...
        %        - conc_Nlab_numerical(t) * (miner + dissol_2);
        %    conc_DIN_numerical(t + 1) = conc_DIN_numerical(t)...
        %        + conc_Nlab_numerical(t) * (miner)...
        %        - conc_DIN_numerical(t) * (denit + plantup)
        %end
        
    end

    %figure
    %plot(c_a, 'k-', 'linewidth',1.5)
    %hold on
    %plot(c_b, 'k--', 'linewidth',1.5)
    %hold on
    %plot(c_c, 'k:', 'linewidth',1.5)
    %xlabel("Time")
    %ylabel("Concentration")
    %legend("Species A",...
    %        "Species B",...
    %        "Species C")
    %grid on
    
    %figure
    %subplot(1,3,1)
    %plot(conc_Nlab_numerical)
    %hold on
    %plot(conc_Nlab_analytical)
    %subplot(1,3,2)
    %plot(conc_Nref_numerical)
    %hold on
    %plot(conc_Nref_analytical)
    %subplot(1,3,3)
    %plot(conc_DIN_numerical)
    %hold on
    %plot(conc_DIN_analytical)
    
end

% Streeter_Phelps -> DO sag
function [c_bod, c_do] = test_Streeter_Phelps(c0_bod, c0_do, k_bod, k_rear, tsim)

    c_bod = zeros(tsim, 1);
    c_do = zeros(tsim, 1);
    
    % Species 1 (simple single species linear decay)
    for t = 1:tsim
        
        ti = t - 1;
        % BOD
        c_bod(t) = c0_bod * exp(-k_bod * ti);
        
        % DO
        A = k_bod * c0_bod / (k_rear - k_bod);
        B = exp(-k_bod * ti) - exp(-k_rear * ti);
        C = c0_do * exp(-k_rear * ti);

        c_do(t) = A * B + C;
    
    end

end

% Plot results (1 species)
function plotAnalytical_singleSpec(conc_A_analytical,...
                                   conc_A_openwq,...
                                   title_str)
   
   time = conc_A_openwq.Time;
   figure
   plot(time, conc_A_analytical,'k-','linewidth',1)
   hold on
   plot(time, conc_A_openwq.data_save_final,'k--s','linewidth',1)
   
    xlabel("Time")
    ylabel("Concentration")
    
    legend("Analytical",...
            "OpenWQ")
        
    title(title_str)
    
    grid on
    
    tsim =  875;
    k_a = 0.01;
    c_a_numerical = zeros(tsim, 1);
    c0_a = 10;
    c_a_numerical(1) = c0_a;
    for t = 1:tsim-1
        
        ti = t - 1;
        
    % Numerical approximation for testing purposes
        dca_dt = k_a * c_a_numerical(t)^2 * (t-ti);
        c_a_numerical(t+1) = c_a_numerical(t) - dca_dt;
    end
     hold on
     plot(time, c_a_numerical,'r-','linewidth',1)
     
end


% Plot results (1 species)
function plotAnalytical_2species(conc_A_analytical,...
                                 conc_B_analytical,...
                                 conc_A_openwq,...
                                 conc_B_openwq,...
                                 title_str,...
                                 spec_name_A,...
                                 spec_name_B...
                                 )
    
    time = conc_A_openwq.Time;
                             
    figure
    sgtitle(title_str)
    
    subplot(1,2,1)
    plot(time, conc_A_analytical,'k-','linewidth',1)
    hold on
    plot(time, conc_A_openwq.data_save_final,'k--s','linewidth',1)
    xlabel("Time")
    ylabel("Concentration")
    legend("Analytical Solution",...
            "OpenWQ")
    title(spec_name_A)
    grid on

    subplot(1,2,2)
    plot(time, conc_B_analytical,'k-','linewidth',1.5)
    hold on
    plot(time, conc_B_openwq.data_save_final,'k--s','linewidth',1)

    xlabel("Time")
    ylabel("Concentration")

    legend("Analytical Solution",...
            "OpenWQ")

    title(spec_name_B)
    grid on
    
end

function plotAnalytical_3species(conc_A_analytical,...
                                 conc_B_analytical,...
                                 conc_C_analytical,...
                                 conc_A_openwq,...
                                 conc_B_openwq,...
                                 conc_C_openwq,...
                                 title_str,...
                                 spec_name_A,...
                                 spec_name_B,...
                                 spec_name_C)
    
   time = conc_A_openwq.Time;
   
    figure
    sgtitle(title_str)
    
    subplot(1,3,1)
    plot(time, conc_A_analytical,'k-','linewidth',1)
    hold on
    plot(time, conc_A_openwq.data_save_final,'k--s','linewidth',1)
    xlabel("Time")
    ylabel("Concentration")
    legend("Analytical Solution",...
            "OpenWQ")
    title(spec_name_A)
    grid on

    subplot(1,3,2)
    plot(time, conc_B_analytical,'k-','linewidth',1.5)
    hold on
    plot(time, conc_B_openwq.data_save_final,'k--s','linewidth',1)
    xlabel("Time")
    ylabel("Concentration")
    legend("Analytical Solution",...
            "OpenWQ")
    title(spec_name_B)
    grid on
        
    subplot(1,3,3)
    plot(time, conc_C_analytical,'k-','linewidth',1.5)
    hold on
    plot(time, conc_C_openwq.data_save_final,'k--s','linewidth',1)
    xlabel("Time")
    ylabel("Concentration")
    legend("Analytical Solution",...
            "OpenWQ")
    title(spec_name_C)
    grid on
    
    
end

function plotAnalytical_NcycleRes(conc_Nref_analytical,...
                            conc_Nlab_analytical,...
                            conc_DON_analytical,...
                            conc_DIN_analytical,...
                            conc_Nref_openwq,...
                            conc_Nlab_openwq,...
                            conc_DON_openwq,...
                            conc_DIN_openwq,...
                            title_str,...
                            Nref_name,...
                            Nlab_name,...
                            DON_name,...
                            DIN_name)
     
                        
    time = conc_Nref_openwq.Time;
   
    figure
    sgtitle(title_str)
    
    subplot(2,2,1)
    plot(time, conc_Nref_analytical,'k-','linewidth',1)
    hold on
    plot(time, conc_Nref_openwq.data_save_final,'k--s','linewidth',1)
    xlabel("Time")
    ylabel("Concentration")
    legend("Analytical Solution",...
            "OpenWQ")
    title(Nref_name)
    grid on
    
    subplot(2,2,2)
    plot(time, conc_Nlab_analytical,'k-','linewidth',1)
    hold on
    plot(time, conc_Nlab_openwq.data_save_final,'k--s','linewidth',1)
    xlabel("Time")
    ylabel("Concentration")
    legend("Analytical Solution",...
            "OpenWQ")
    title(Nlab_name)
    grid on
    
    subplot(2,2,3)
    plot(time, conc_DON_analytical,'k-','linewidth',1)
    hold on
    plot(time, conc_DON_openwq.data_save_final,'k--s','linewidth',1)
    xlabel("Time")
    ylabel("Concentration")
    legend("Analytical Solution",...
            "OpenWQ")
    title(DON_name)
    grid on
    
    subplot(2,2,4)
    plot(time, conc_DIN_analytical,'k-','linewidth',1)
    hold on
    plot(time, conc_DIN_openwq.data_save_final,'k--s','linewidth',1)
    xlabel("Time")
    ylabel("Concentration")
    legend("Analytical Solution",...
            "OpenWQ")
    title(DIN_name)
    grid on
    
end