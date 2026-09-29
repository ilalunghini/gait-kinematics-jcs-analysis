%% 01 | SETUP AND DATA LOADING
% Load marker trajectories and gait-cycle events from acquisition.mat.
        clear; close all; clc;
        fprintf('\nStarting gait kinematics analysis...\n');


        scriptDirectory = fileparts(mfilename('fullpath'));
        repositoryDirectory = fileparts(scriptDirectory);
        dataFile = fullfile(repositoryDirectory, 'data', 'acquisition.mat');

        if ~isfile(dataFile)
            try
                activeFile = matlab.desktop.editor.getActiveFilename;
                activeScriptDirectory = fileparts(activeFile);
                activeRepositoryDirectory = fileparts(activeScriptDirectory);
                activeDataFile = fullfile(activeRepositoryDirectory, 'data', 'acquisition.mat');
                if isfile(activeDataFile)
                    scriptDirectory = activeScriptDirectory;
                    repositoryDirectory = activeRepositoryDirectory;
                    dataFile = activeDataFile;
                end
            catch
            end
        end

        addpath(fullfile(scriptDirectory, 'functions'));

        if ~isfile(dataFile)
            error('Required data file not found: %s', dataFile);
        end

        fprintf('[1/7] Loading marker trajectories and gait-cycle events...\n');
        load(dataFile, 'static', 'dynamic', 'GaitCycle');
        fprintf('      Loaded %d static frames, %d dynamic frames, and %d gait cycles.\n', ...
            size(static.rp5m,1), size(dynamic.rmma,1), size(GaitCycle,1));
        % Static marker trajectories.
        RMMA = static.rmma;
        RANK = static.rank;
        RHFB = static.rhfb;
        RCA1 = static.rca1;
        RCA2 = static.rca2;
        RLCA = static.rlca;
        RSTL = static.rstl;
        RD1M = static.rd1m;
        RD5M = static.rd5m;
        RP5M = static.rp5m;
        RP1M = static.rp1m;
        RTOE = static.rtoe;
        RTUB = static.rtub;

        nF_st = size(RP5M,1);
       %% 02 | STATIC CALIBRATION AND LOCAL COORDINATE SYSTEMS
        fprintf('[2/7] Building static segment frames and marker templates...\n');
       % Tibia.
        O_tibia = (RMMA + RANK) / 2;
        z_tibia = (RANK - RMMA) ./ vecnorm(RANK - RMMA,2,2); 
        x_tibia = cross(RHFB - O_tibia, z_tibia, 2);
        x_tibia = x_tibia ./ vecnorm(x_tibia,2,2);
        y_tibia = cross(z_tibia,x_tibia);
        y_tibia = y_tibia ./ vecnorm(y_tibia,2,2);
        T_tibia = repmat(eye(4),[1, 1, nF_st]);
        T_tibia(1:3,1,:) = permute(x_tibia, [2, 3, 1]);
        T_tibia(1:3,2,:) = permute(y_tibia, [2, 3, 1]);
        T_tibia(1:3,3,:) = permute(z_tibia, [2, 3, 1]);
        T_tibia(1:3,4,:) = permute(O_tibia, [2, 3, 1]);
        
        RTUB_tb = mean(coordChange(Tinv(T_tibia), static.rtub), 1);
        RANK_tb = mean(coordChange(Tinv(T_tibia), static.rank), 1);
        RHFB_tb = mean(coordChange(Tinv(T_tibia), static.rhfb), 1);
        RMMA_tb = mean(coordChange(Tinv(T_tibia), static.rmma), 1);
        RSHN_tb = mean(coordChange(Tinv(T_tibia), static.rshn), 1);
        
        Points_tb = [RTUB_tb; RANK_tb; RHFB_tb; RMMA_tb; RSHN_tb];
        Points_tb = Points_tb - mean(Points_tb, 1);

        % Calcaneus.
        O_calcagno = (RCA1 + RCA2 + RLCA + RSTL) / 4; 
        punto_medio = (RLCA + RSTL) / 2;
        normale_piano = cross(RCA2 - RCA1, RCA1 - punto_medio, 2); 
        z_calcagno = normale_piano ./ vecnorm(normale_piano,2,2); 
        z_piano_appoggio = repmat([0,0,1], [nF_st, 1]); 
        x_calcagno = cross(z_piano_appoggio, z_calcagno,2);  
        x_calcagno = x_calcagno ./ vecnorm(x_calcagno,2,2);  
        y_calcagno = cross(z_calcagno, x_calcagno,2);
        y_calcagno = y_calcagno ./ vecnorm(y_calcagno,2,2);
        
        T_calcagno = repmat(eye(4),[1, 1, nF_st]);
        T_calcagno(1:3,1,:) = permute(x_calcagno, [2, 3, 1]);
        T_calcagno(1:3,2,:) = permute(y_calcagno, [2, 3, 1]);
        T_calcagno(1:3,3,:) = permute(z_calcagno, [2, 3, 1]);
        T_calcagno(1:3,4,:) = permute(O_calcagno, [2, 3, 1]); 

        RCA1_ca = mean(coordChange(Tinv(T_calcagno), static.rca1), 1);
        RCA2_ca = mean(coordChange(Tinv(T_calcagno), static.rca2), 1);
        RLCA_ca = mean(coordChange(Tinv(T_calcagno), static.rlca), 1);
        RSTL_ca = mean(coordChange(Tinv(T_calcagno), static.rstl), 1);
        
        Points_ca = [RCA1_ca; RCA2_ca; RLCA_ca; RSTL_ca];
        Points_ca = Points_ca - mean(Points_ca, 1);

        % Metatarsus.
        O_metatarso = (RD1M + RD5M + RP5M) / 3;
        y_metatarso = cross(RD5M - RP5M, RD1M - RP5M, 2);
        y_metatarso = y_metatarso ./ vecnorm(y_metatarso, 2, 2);
        x_metatarso = RTOE - ((RP1M + RP5M) / 2); 
        x_metatarso = x_metatarso - (dot(x_metatarso, y_metatarso, 2) .* y_metatarso);
        x_metatarso = x_metatarso ./ vecnorm(x_metatarso, 2, 2);
        z_metatarso = cross(x_metatarso, y_metatarso, 2);

        
        T_metatarso = repmat(eye(4),[1, 1, nF_st]);
        T_metatarso(1:3,1,:) = permute(x_metatarso, [2, 3, 1]);
        T_metatarso(1:3,2,:) = permute(y_metatarso, [2, 3, 1]);
        T_metatarso(1:3,3,:) = permute(z_metatarso, [2, 3, 1]);
        T_metatarso(1:3,4,:) = permute(O_metatarso, [2, 3, 1]); 

        RP1M_me = mean(coordChange(Tinv(T_metatarso), static.rp1m), 1);
        RP5M_me = mean(coordChange(Tinv(T_metatarso), static.rp5m), 1);
        RD1M_me = mean(coordChange(Tinv(T_metatarso), static.rd1m), 1);
        RD5M_me = mean(coordChange(Tinv(T_metatarso), static.rd5m), 1);
        RTOE_me = mean(coordChange(Tinv(T_metatarso), static.rtoe), 1);
        
        % Exclude RD1M_me because RD1M is unavailable in the dynamic trial.
        Points_me = [RP5M_me; RP1M_me; RD5M_me; RTOE_me];
        Points_me = Points_me - mean(Points_me, 1);
        

        % Technical foot.
        O_piede_tecnico = RCA2;
        % The assignment defines technical X from RP5M to RP1M.
        x_piede_tecnico = (RP1M - RP5M) ./ vecnorm(RP1M - RP5M, 2, 2);
        y_piede_tecnico = cross(x_piede_tecnico, RCA2 - RP5M, 2);
        y_piede_tecnico= y_piede_tecnico ./ vecnorm(y_piede_tecnico,2,2);
        z_piede_tecnico = cross(x_piede_tecnico, y_piede_tecnico, 2) ./ vecnorm(cross(x_piede_tecnico, y_piede_tecnico, 2), 2, 2); 
       

        T_piede_tecnico = repmat(eye(4),[1, 1, nF_st]);
        T_piede_tecnico(1:3,1,:) = permute(x_piede_tecnico, [2, 3, 1]);
        T_piede_tecnico(1:3,2,:) = permute(y_piede_tecnico, [2, 3, 1]);
        T_piede_tecnico(1:3,3,:) = permute(z_piede_tecnico, [2, 3, 1]);
        T_piede_tecnico(1:3,4,:) = permute(O_piede_tecnico, [2, 3, 1]); 

        % These local-marker trajectories are used to construct
        % the dynamic reference frames.
        RCA2_piede = coordChange(Tinv(T_piede_tecnico), static.rca2);  
        RCA1_piede = coordChange(Tinv(T_piede_tecnico), static.rca1);  
        RP1M_piede = coordChange(Tinv(T_piede_tecnico), static.rp1m);  
        RP5M_piede = coordChange(Tinv(T_piede_tecnico), static.rp5m);  
        RTOE_piede = coordChange(Tinv(T_piede_tecnico), static.rtoe);  
        Points_piede_tc = [RCA1_piede; RP1M_piede; RP5M_piede; RTOE_piede];
        Points_piede_tc = Points_piede_tc - mean(Points_piede_tc, 1);
        


        % Anatomical foot defined relative to the technical foot.
        O_piede_anatomico = RCA2_piede;  
        x_piede_anatomico = (RTOE_piede - RCA2_piede) ./ vecnorm(RTOE_piede - RCA2_piede, 2, 2);
        % All vectors here are expressed in technical-foot coordinates.
        % Technical Y is normal to the technical XZ plane.
        technicalY = repmat([0, 1, 0], nF_st, 1);
        technicalX = repmat([1, 0, 0], nF_st, 1);
        x_piede_anatomico = x_piede_anatomico - dot(x_piede_anatomico, technicalY, 2) .* technicalY;
        x_piede_anatomico = x_piede_anatomico ./ vecnorm(x_piede_anatomico, 2, 2);  
        y_piede_anatomico = cross(x_piede_anatomico, technicalX, 2);
        y_piede_anatomico = y_piede_anatomico ./ vecnorm(y_piede_anatomico, 2, 2);  
        z_piede_anatomico = cross(y_piede_anatomico, x_piede_anatomico, 2);  
        z_piede_anatomico = -z_piede_anatomico ./ vecnorm(z_piede_anatomico, 2, 2);

        % Anatomical-foot transformation relative to the technical foot.
        T_tec_piede_anatomico = repmat(eye(4),[1, 1, nF_st]);
        T_tec_piede_anatomico(1:3,1,:) = permute(x_piede_anatomico, [2, 3, 1]);
        T_tec_piede_anatomico(1:3,2,:) = permute(y_piede_anatomico, [2, 3, 1]);
        T_tec_piede_anatomico(1:3,3,:) = permute(z_piede_anatomico, [2, 3, 1]);
        T_tec_piede_anatomico(1:3,4,:) = permute(O_piede_anatomico, [2, 3, 1]); 

        %% Transform the anatomical-foot frame from technical-foot to laboratory coordinates.
        T_piede_anatomico = pagemtimes(T_piede_tecnico, T_tec_piede_anatomico);
    
        %% Average foot-marker coordinates for dynamic registration.
        RCA2_piede= mean(coordChange(Tinv(T_piede_anatomico), static.rca2) , 1);
        RCA1_piede = mean(coordChange(Tinv(T_piede_anatomico), static.rca1), 1);
        RP1M_piede = mean(coordChange(Tinv(T_piede_anatomico), static.rp1m), 1);
        RP5M_piede = mean(coordChange(Tinv(T_piede_anatomico), static.rp5m), 1);
        RTOE_piede = mean(coordChange(Tinv(T_piede_anatomico), static.rtoe), 1);
        
        Points_piede_an = [RCA1_piede; RP1M_piede; RP5M_piede; RTOE_piede];
        Points_piede_an = Points_piede_an - mean(Points_piede_an, 1);

       
        %% 03 | STATIC FRAME VISUALIZATION
        fprintf('[3/7] Plotting static calibration frames...\n');
        %figure;
        hold on;

        frame = randi([1, 641]);

         plotCS(T_tibia(:, :, frame), 200);
         plotCS(T_calcagno(:, :, frame), 200);
         plotCS(T_metatarso(:, :, frame), 200); 
         plotCS(T_piede_anatomico(:, :, frame), 200);

        axis equal;
        
        % Plot static tibia markers.
        % RHFB
        scatter3(static.rhfb(:,1), static.rhfb(:,2), static.rhfb(:,3), 'k', 'DisplayName', 'RHFB');
        text(static.rhfb(1,1), static.rhfb(1,2), static.rhfb(1,3),  'RHFB', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal');

        % RTUB
        scatter3(static.rtub(:,1), static.rtub(:,2), static.rtub(:,3), 'k', 'DisplayName', 'RTUB');
        text(static.rtub(1,1), static.rtub(1,2), static.rtub(1,3), 'RTUB', 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');

        % RMMA
        scatter3(static.rmma(:,1), static.rmma(:,2), static.rmma(:,3), 'k', 'DisplayName', 'RMMA');
        text(static.rmma(1,1), static.rmma(1,2), static.rmma(1,3), 'RMMA', 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');

        % RANK
        scatter3(static.rank(:,1), static.rank(:,2), static.rank(:,3), 'k', 'DisplayName', 'RANK');
        text(static.rank(1,1), static.rank(1,2), static.rank(1,3), 'RANK', 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');

        % scatter3(static.rshn(:,1), static.rshn(:,2), static.rshn(:,3),
        % 'ro', 'DisplayName', 'RSHN'); omit this marker to keep the
        % illustrative polygon regular.

        % Polygon vertices at the selected frame.
        p1 = static.rhfb(frame, :); 
        p2 = static.rtub(frame, :); 
        p3 = static.rmma(frame, :); 
        p4 = static.rank(frame, :); 
        %p5 = static.rshn(frame, :); 

        % Polygon-point matrix.
        points = [p1; p2; p3; p4; p1]; 

        fill3(points(:,1), points(:,2), points(:,3), 'k', 'FaceAlpha', 0.2); 

        % Plot static calcaneus markers.
        % RCA1
        scatter3(static.rca1(:,1), static.rca1(:,2), static.rca1(:,3), 'k', 'DisplayName', 'RCA1');
        text(static.rca1(1,1), static.rca1(1,2), static.rca1(1,3) , 'RCA1', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % Omit RCA2 to keep the illustrative polygon regular.
        % scatter3(static.rca2(:,1), static.rca2(:,2), static.rca2(:,3), 'k', 'DisplayName', 'RCA2');
        % text(static.rca2(1,1), static.rca2(1,2), static.rca2(1,3) + 5, 'RCA2', ...
        % 'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % RLCA
        scatter3(static.rlca(:,1), static.rlca(:,2), static.rlca(:,3), 'k', 'DisplayName', 'RLCA');
        text(static.rlca(1,1), static.rlca(1,2), static.rlca(1,3) , 'RLCA', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % RSTL
        scatter3(static.rstl(:,1), static.rstl(:,2), static.rstl(:,3), 'k', 'DisplayName', 'RSTL');
        text(static.rstl(1,1), static.rstl(1,2), static.rstl(1,3) , 'RSTL', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % Polygon vertices at the selected frame.
        p1 = static.rca1(frame, :);
        % p2 = static.rca2(frame, :);
        p3 = static.rstl(frame, :);
        p4 = static.rlca(frame, :);

        % Polygon-point matrix.
        points = [p1; p3; p4; p1];

        % Draw a transparent surface.
        fill3(points(:,1), points(:,2), points(:,3), 'g', 'FaceAlpha', 0.2);

        % Plot static metatarsus markers.
        % RP1M
        scatter3(static.rp1m(:,1), static.rp1m(:,2), static.rp1m(:,3), 'k', 'DisplayName', 'RP1M');
        text(static.rp1m(1,1), static.rp1m(1,2), static.rp1m(1,3) + 5, 'RP1M', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % RP5M
        scatter3(static.rp5m(:,1), static.rp5m(:,2), static.rp5m(:,3), 'k', 'DisplayName', 'RP5M');
        text(static.rp5m(1,1), static.rp5m(1,2), static.rp5m(1,3) + 5, 'RP5M', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % Omit RD1M to keep the illustrative polygon regular.
        % scatter3(static.rd1m(:,1), static.rd1m(:,2), static.rd1m(:,3), 'k', 'DisplayName', 'RD1M');
        % text(static.rd1m(1,1), static.rd1m(1,2), static.rd1m(1,3) + 5, 'RD1M', ...
        %     'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');
        % 
        % Omit RD5M to keep the illustrative polygon regular.
        % scatter3(static.rd5m(:,1), static.rd5m(:,2), static.rd5m(:,3), 'k', 'DisplayName', 'RD5M');
        % text(static.rd5m(1,1), static.rd5m(1,2), static.rd5m(1,3) + 5, 'RD5M', ...
        %     'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');
        % 

        % RTOE
        scatter3(static.rtoe(:,1), static.rtoe(:,2), static.rtoe(:,3), 'k', 'DisplayName', 'RTOE');
        text(static.rtoe(1,1), static.rtoe(1,2), static.rtoe(1,3) + 5, 'RTOE', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % Polygon vertices at the selected frame.
        p1 = static.rp1m(frame, :);  
        p2 = static.rp5m(frame, :);  
        %p3 = static.rd5m(frame, :);  
        %p4 = static.rd1m(frame, :);  
        p5 = static.rtoe(frame, :);

        % Polygon-point matrix.
        points = [p1; p2; p5; p1];

        fill3(points(:,1), points(:,2), points(:,3), 'r', 'FaceAlpha', 0.2);

        
        % Plot static anatomical-foot markers.
        % RCA1
        scatter3(static.rca1(:,1), static.rca1(:,2), static.rca1(:,3), 'k', 'DisplayName', 'RCA1');
        text(static.rca1(1,1), static.rca1(1,2), static.rca1(1,3) + 5, 'RCA1', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % RCA2
        scatter3(static.rca2(:,1), static.rca2(:,2), static.rca2(:,3), 'k', 'DisplayName', 'RCA1');
        text(static.rca2(1,1), static.rca2(1,2), static.rca2(1,3) + 5, 'RCA2', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % RP1M
        scatter3(static.rp1m(:,1), static.rp1m(:,2), static.rp1m(:,3), 'k', 'DisplayName', 'RP1M');
        text(static.rp1m(1,1), static.rp1m(1,2), static.rp1m(1,3) + 5, 'RP1M', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % RP5M
        scatter3(static.rp5m(:,1), static.rp5m(:,2), static.rp5m(:,3), 'k', 'DisplayName', 'RP5M');
        text(static.rp5m(1,1), static.rp5m(1,2), static.rp5m(1,3) + 5, 'RP5M', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');

        % RTOE
        scatter3(static.rtoe(:,1), static.rtoe(:,2), static.rtoe(:,3), 'k', 'DisplayName', 'RTOE');
        text(static.rtoe(1,1), static.rtoe(1,2), static.rtoe(1,3) + 5, 'RTOE', ...
            'FontSize', 10, 'Color', 'k', 'FontWeight', 'normal', 'VerticalAlignment', 'bottom');



        % Polygon vertices at the selected frame.
        %p1 = static.rca1(frame, :); % Omit to keep the illustrative polygon regular.
        p2 = static.rca2(frame, :);
        p3 = static.rp1m(frame, :);   
        p4 = static.rtoe(frame, :);   
        p5 = static.rp5m(frame, :);   

        % Polygon-point matrix.
        points = [ p2; p3; p4; p5; p2];

        fill3(points(:,1), points(:,2), points(:,3), 'b', 'FaceAlpha', 0.2);

        xlabel('X (red)'); ylabel('Y (green)'); zlabel('Z (blue)');
        title('Static Local Coordinate Systems');
        grid on;
        view(3);
        hold off;
        
        % Average technical-foot marker coordinates for dynamic registration.
        RCA2_piede = mean(coordChange(Tinv(T_piede_tecnico), static.rca2));  
        RCA1_piede = mean(coordChange(Tinv(T_piede_tecnico), static.rca1));  
        RP1M_piede = mean(coordChange(Tinv(T_piede_tecnico), static.rp1m));  
        RP5M_piede = mean(coordChange(Tinv(T_piede_tecnico), static.rp5m));  
        RTOE_piede = mean(coordChange(Tinv(T_piede_tecnico), static.rtoe));  
        Points_piede_tc = [RCA1_piede; RP1M_piede; RP5M_piede; RTOE_piede];
        Points_piede_tc = Points_piede_tc - mean(Points_piede_tc, 1);

        %% 04 | DYNAMIC SEGMENT REGISTRATION
        fprintf('[4/7] Estimating dynamic segment rotations using SVD...\n');
        % Total number of frames, obtained from a marker trajectory in the
        % dynamic structure.
        numFrames = size(dynamic.rmma, 1);

        % Points_tb = [RTUB_tb; RANK_tb; RHFB_tb; RMMA_tb; RSHN_tb];
        Points_tb_lab = permute(cat(3, dynamic.rtub, dynamic.rank, dynamic.rhfb, dynamic.rmma, dynamic.rshn), [2, 3, 1]);
        % Center markers at each frame.
        Points_tb_lab = Points_tb_lab - mean(Points_tb_lab, 2, "omitmissing");

        % Points_ca = [RCA1_ca; RCA2_ca; RLCA_ca; RSTL_ca];
        Points_ca_lab=permute(cat(3, dynamic.rca1, dynamic.rca2, dynamic.rlca, dynamic.rstl), [2,3,1]);
        % Center markers at each frame.
        Points_ca_lab = Points_ca_lab - mean(Points_ca_lab, 2, "omitmissing");

        % Points_me = [RP1M_me; RP5M_me; RD5M_me; RTOE_me]; 
        % RD1M is unavailable in the dynamic trial. It is excluded from both
        % static and dynamic registration while preserving marker order.

        Points_me_lab=permute(cat(3, dynamic.rp5m, dynamic.rp1m, dynamic.rd5m, dynamic.rtoe), [2,3,1]); 
        % Center markers at each frame.
        Points_me_lab=Points_me_lab - mean(Points_me_lab,2, "omitmissing");


        % Points_piede_tc (and anatomical) = [RCA1_piede; RP1M_piede; RP5M_piede;
        % RTOE_piede]; 
        Points_piede_tc_lab = permute(cat(3, dynamic.rca1, dynamic.rp1m, dynamic.rp5m, dynamic.rtoe), [2,3,1]);
        % Center technical-foot markers at each frame.
        Points_piede_tc_lab =  Points_piede_tc_lab - mean(Points_piede_tc_lab, 2, "omitmissing");
        Points_piede_an_lab = permute(cat(3, dynamic.rca1, dynamic.rp1m, dynamic.rp5m, dynamic.rtoe), [2,3,1]);
        % Center anatomical-foot markers at each frame.
        Points_piede_an_lab = Points_piede_an_lab - mean(Points_piede_an_lab, 2 , "omitmissing");

        % Compute cross-covariance matrices for each segment.
        Cross_tb = CrossProduct(Points_tb_lab, Points_tb);
        Cross_ca = CrossProduct(Points_ca_lab, Points_ca);
        Cross_me = CrossProduct(Points_me_lab, Points_me);
        Cross_piede_tc = CrossProduct(Points_piede_tc_lab, Points_piede_tc);
        Cross_piede_an = CrossProduct(Points_piede_an_lab, Points_piede_an);

        % Compute least-squares rotation matrices for each segment.
        R_tb = computeRotationMatrix(Cross_tb);
        R_ca = computeRotationMatrix(Cross_ca);
        R_me = computeRotationMatrix(Cross_me);
        R_piede_tc = computeRotationMatrix(Cross_piede_tc);
        R_piede_an = computeRotationMatrix(Cross_piede_an);

        % Compute dynamic segment origins.
        origin_tb       = permute((dynamic.rmma + dynamic.rank) / 2, [2, 3, 1]);    
        origin_ca       = permute((dynamic.rca1 + dynamic.rca2 + dynamic.rlca + dynamic.rstl)/4, [2, 3, 1]);
        origin_me       = permute((dynamic.rd5m + dynamic.rp5m + dynamic.rtoe + dynamic.rp1m) / 4, [2, 3, 1]); % RD1M is unavailable in the dynamic trial.
        origin_piede_tc = permute(dynamic.rca2, [2, 3, 1]);
        origin_piede_an = permute(dynamic.rca2, [2, 3, 1]);


        % Preallocate homogeneous transformations.
        T_tb         = repmat(eye(4), [1, 1, numFrames]);
        T_ca         = repmat(eye(4), [1, 1, numFrames]);
        T_me         = repmat(eye(4), [1, 1, numFrames]);
        T_piede_tc   = repmat(eye(4), [1, 1, numFrames]);
        T_piede_an   = repmat(eye(4), [1, 1, numFrames]);
        
        % Insert rotation matrices.
        T_tb(1:3,1:3,:)         = R_tb;
        T_ca(1:3,1:3,:)         = R_ca;
        T_me(1:3,1:3,:)         = R_me;
        T_piede_tc(1:3,1:3,:)   = R_piede_tc;
        T_piede_an(1:3,1:3,:)   = R_piede_an;
        
        % Insert translations (the dynamic segment origins).
        T_tb(1:3,4,:) = origin_tb;  % Convert [numFrames,3] to [3,1,numFrames].
        T_ca(1:3,4,:)         = origin_ca;
        T_me(1:3,4,:)         = origin_me;
        T_piede_tc(1:3,4,:)   = origin_piede_tc;
        T_piede_an(1:3,4,:)   = origin_piede_an;
        
        %% 05 | DYNAMIC POSE VISUALIZATION
        fprintf('[5/7] Plotting reconstructed dynamic poses...\n');
        framesToPlot = round(linspace(1, numFrames, 4));
        nFrames = length(framesToPlot);
        
        T = {T_tb, T_ca, T_me, T_piede_an};
        segmentNames = {'Tibia', 'Calcaneus', 'Metatarsus', 'Anatomical Foot'};
        nSegments = length(T);
       
        figure('Name', 'Dynamic Coordinate Systems - Separate and Combined', ...
               'Units', 'normalized', 'Position', [0.05 0.4 0.9 0.5]);
        
        % Plot individual segments.
        for s = 1:nSegments
            subplot(1, nSegments + 1, s);  % The final panel is the combined plot.
            hold on;
            axis equal;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            grid on;
            view(3);
            title(segmentNames{s});
        
            for i = 1:nFrames
                f = framesToPlot(i);
                plotCS(T{s}(:,:,f), 400);
            end
        end
        
        % Plot all segments together.
        subplot(1, nSegments + 1, nSegments + 1);
        hold on;
        axis equal;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        grid on;
        view(3);
        title('Coordinate Systems of All Segments');
        
        for i = 1:nFrames
            f = framesToPlot(i);
            plotCS(T_tb(:,:,f), 400);
            plotCS(T_ca(:,:,f), 400);
            plotCS(T_me(:,:,f), 400);
            plotCS(T_piede_tc(:,:,f), 400);
            plotCS(T_piede_an(:,:,f), 400);
        end
        
        hold off;


            %% 06 | JOINT ANGLES AND GAIT-CYCLE NORMALIZATION
            fprintf('[6/7] Calculating joint angles and normalizing gait cycles...\n');
            % Relative rotations map distal-segment coordinates into the
            % proximal-segment coordinate system.
            R_tibia_calcagno = pagemtimes(R_tb,'transpose', R_ca,'none');
            R_calcagno_metatarso = pagemtimes(R_ca,'transpose' , R_me,'none');
            R_tibia_piede = pagemtimes(R_tb,'transpose' ,R_piede_an, 'none');

            % The JCS uses the proximal Z-axis (e1), the distal Y-axis (e3),
            % and their common perpendicular floating axis. The following
            % Z-X-Y decomposition gives alpha about e1, beta about the
            % floating axis, and gamma about e3. No division by cos(beta) is
            % part of these angle definitions.
            alpha_tc = squeeze(atan2(-R_tibia_calcagno(1,2,:), R_tibia_calcagno(2,2,:)));
            beta_tc  = squeeze(asin(R_tibia_calcagno(3,2,:)));
            gamma_tc = squeeze(atan2(-R_tibia_calcagno(3,1,:), R_tibia_calcagno(3,3,:)));

            alpha_cm = squeeze(atan2(-R_calcagno_metatarso(1,2,:), R_calcagno_metatarso(2,2,:)));
            beta_cm  = squeeze(asin(R_calcagno_metatarso(3,2,:)));
            gamma_cm = squeeze(atan2(-R_calcagno_metatarso(3,1,:), R_calcagno_metatarso(3,3,:)));

            alpha_tp = squeeze(atan2(-R_tibia_piede(1,2,:), R_tibia_piede(2,2,:)));
            beta_tp  = squeeze(asin(R_tibia_piede(3,2,:)));
            gamma_tp = squeeze(atan2(-R_tibia_piede(3,1,:), R_tibia_piede(3,3,:)));

            % The assignment labels the three requested angles as e1, e2,
            % and e3: plantar/dorsiflexion, internal/external rotation, and
            % inversion/eversion, respectively.
            Angoli_rotazione_tibiacalcagno_all = [alpha_tc, gamma_tc, beta_tc];
            Angoli_rotazione_calcagnometatarso_all = [alpha_cm, gamma_cm, beta_cm];
            Angoli_rotazione_tibiapiede_all = [alpha_tp, gamma_tp, beta_tp];

            %% Convert angles to degrees
            Angoli_rotazione_tibiacalcagno_all = rad2deg(Angoli_rotazione_tibiacalcagno_all);
            Angoli_rotazione_calcagnometatarso_all = rad2deg(Angoli_rotazione_calcagnometatarso_all);
            Angoli_rotazione_tibiapiede_all = rad2deg(Angoli_rotazione_tibiapiede_all);
            
               
           %% Tibia-calcaneus
            Angoli_rotazione_tibiacalcagno_all = squeeze(Angoli_rotazione_tibiacalcagno_all);

            % Normalize angles to each gait cycle.
            numCicli = size(GaitCycle, 1);
            angoli_tibiacalcagno_norm = zeros(101, 3, numCicli); % Preallocate normalized curves.
            
            for icycle = 1:numCicli
                angoli_tibiacalcagno_norm(:,:,icycle) = eventsnormalize(Angoli_rotazione_tibiacalcagno_all, GaitCycle(icycle, :), 101);
            end

           %% Calcaneus-metatarsus
            Angoli_rotazione_calcagnometatarso_all = squeeze(Angoli_rotazione_calcagnometatarso_all);

            % Normalize angles to each gait cycle.
            numCicli = size(GaitCycle, 1);
            angoli_calcagnometatarso_norm = zeros(101, 3, numCicli); % Preallocate normalized curves.
            
            for icycle = 1:numCicli
                angoli_calcagnometatarso_norm(:,:,icycle) = eventsnormalize(Angoli_rotazione_calcagnometatarso_all, GaitCycle(icycle, :), 101);
            end

             
           %% Tibia-anatomical foot (Ankle 2)
            Angoli_rotazione_tibiapiede_all = squeeze(Angoli_rotazione_tibiapiede_all);

            % Normalize angles to each gait cycle.
            numCicli = size(GaitCycle, 1);
            angoli_tibiapiede_norm = zeros(101, 3, numCicli); % Preallocate normalized curves.
            
            for icycle = 1:numCicli
                angoli_tibiapiede_norm(:,:,icycle) = eventsnormalize(Angoli_rotazione_tibiapiede_all, GaitCycle(icycle, :), 101);
            end

           
            normalizedAngles = {angoli_tibiacalcagno_norm, ...
                                angoli_calcagnometatarso_norm, ...
                                angoli_tibiapiede_norm};
            jointNames = {'Ankle 1: Tibia-Calcaneus', ...
                          'Calcaneus-Metatarsus', ...
                          'Ankle 2: Tibia-Anatomical Foot'};
            angleNames = {'Plantar/Dorsiflexion (e1)', ...
                          'Internal/External Rotation (e2)', ...
                          'Inversion/Eversion (e3)'};
            gaitPercent = linspace(0, 100, 101)';

            % Plot all normalized joint angles. Each panel contains one line
            % per detected gait cycle and is fully labelled as required.
            figure('Name', 'Normalized Joint Angles');
            for jointIndex = 1:3
                for angleIndex = 1:3
                    subplot(3,3,(jointIndex - 1) * 3 + angleIndex); hold on;
                    plot(gaitPercent, squeeze(normalizedAngles{jointIndex}(:, angleIndex, :)));
                    title([jointNames{jointIndex} ' - ' angleNames{angleIndex}]);
                    xlabel('Gait cycle (%)');
                    ylabel('Angle (deg)');
                    grid on;
                end
            end

 
 
                
                %% 07 | SAGITTAL ANGULAR VELOCITY ANALYSIS
                fprintf('[7/7] Calculating and plotting sagittal angular velocities...\n');
                % Sampling frequency in Hz. This value must match the data
                % acquisition setting.
                freq = 100;
                fprintf('      Using the configured sampling frequency: %g Hz.\n', freq);
                dt = 1 / freq;           
                
                % Number of dynamic frames.
                nFrames = size(R_tb, 3); 
                
                % Preallocate three-component angular-velocity vectors.
                omega_tb = zeros(3, nFrames);
                omega_foot = zeros(3, nFrames);
                omega_ankle2 = zeros(3, nFrames);
                
                
                for i = 2:nFrames-1
                    % Central numerical derivative of each rotation matrix.
                    dR_tb = (R_tb(:,:,i+1) - R_tb(:,:,i-1)) / (2 * dt);
                    dR_foot = (R_piede_an(:,:,i+1) - R_piede_an(:,:,i-1)) / (2 * dt);
                    dR_ankle2 = (R_tibia_piede(:,:,i+1) - R_tibia_piede(:,:,i-1)) / (2 * dt);
                
                    % Skew-symmetric angular-velocity matrices. Segment
                    % velocities are expressed in the laboratory frame;
                    % Ankle 2 velocity is expressed in the tibia frame.
                    Omega_tb = dR_tb * R_tb(:,:,i)';
                    Omega_foot = dR_foot * R_piede_an(:,:,i)';
                    Omega_ankle2 = dR_ankle2 * R_tibia_piede(:,:,i)';
                
                    % Angular-velocity vectors.
                    omega_tb(:,i) = [Omega_tb(3,2); Omega_tb(1,3); Omega_tb(2,1)];
                    omega_foot(:,i) = [Omega_foot(3,2); Omega_foot(1,3); Omega_foot(2,1)];
                    omega_ankle2(:,i) = [Omega_ankle2(3,2); Omega_ankle2(1,3); Omega_ankle2(2,1)];
                end
                
                % Convert from radians to degrees per second.
                omega_tb = rad2deg(omega_tb);
                omega_foot = rad2deg(omega_foot);
                omega_ankle2 = rad2deg(omega_ankle2);
                

                % The sagittal component is rotation about the local
                % medio-lateral Z-axis. The segment axes are converted to
                % laboratory coordinates; the Ankle 2 relative velocity is
                % already expressed in tibial coordinates.
                z_tb = squeeze(R_tb(:,3,:));
                z_foot = squeeze(R_piede_an(:,3,:));
                omega_tb_sag = zeros(nFrames,1);
                omega_foot_sag = zeros(nFrames,1);
                for i = 1:nFrames
                    omega_tb_sag(i) = dot(omega_tb(:,i), z_tb(:,i));
                    omega_foot_sag(i) = dot(omega_foot(:,i), z_foot(:,i));
                end
                omega_ankle2_sag = omega_ankle2(3,:)';
                
                %% Normalize angular velocities to the gait cycle
                omega_tb_sag_norm = eventsnormalize(omega_tb_sag, GaitCycle, 101);
                omega_foot_sag_norm = eventsnormalize(omega_foot_sag, GaitCycle, 101);
                omega_ankle2_sag_norm = eventsnormalize(omega_ankle2_sag, GaitCycle, 101);

                % Arrange the concatenated output as one 101-point column
                % per gait cycle before plotting.
                omega_tb_sag_norm = reshape(omega_tb_sag_norm, 101, numCicli);
                omega_foot_sag_norm = reshape(omega_foot_sag_norm, 101, numCicli);
                omega_ankle2_sag_norm = reshape(omega_ankle2_sag_norm, 101, numCicli);
                
                %% Plot normalized sagittal angular velocities
                figure('Name', 'Normalized Sagittal Angular Velocities');

                % Tibia.
                subplot(3,1,1);
                plot(gaitPercent, omega_tb_sag_norm, 'r');
                title('Tibia Sagittal Angular Velocity');
                xlabel('Gait cycle (%)');
                ylabel('Angular velocity (deg/s)');
                grid on;

                % Anatomical foot.
                subplot(3,1,2);
                plot(gaitPercent, omega_foot_sag_norm, 'b');
                title('Anatomical Foot Sagittal Angular Velocity');
                xlabel('Gait cycle (%)');
                ylabel('Angular velocity (deg/s)');
                grid on;

                % Ankle 2: anatomical foot relative to tibia.
                subplot(3,1,3);
                plot(gaitPercent, omega_ankle2_sag_norm, 'g');
                title('Ankle 2 Sagittal Angular Velocity');
                xlabel('Gait cycle (%)');
                ylabel('Angular velocity (deg/s)');
                grid on;
                
 hold off;
 fprintf('\nGait kinematics analysis completed successfully.\n');
