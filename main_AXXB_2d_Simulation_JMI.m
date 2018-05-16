
% /***************************************************************************
% Copyright 
% MUSiiC Laboratory
% Haichong Zhang,Emad M Boctor
% Johns Hopkins University
% 
% For commercial use/licensing, please contact Emad Boctor, Ph.D. at eboctor@jhmi.edu.
% ***************************************************************************/

%% Parameter Setting
% initialization
clear all;
close all;

% pose to be used to 
% eposes = 60;
eposes = 12;

% number of points on phantom
npoint = 12;

% magnitude of segmentation error
MagErrorSeg = 0;

% magnitude out of plane error
MagErrorOOP = 1;

% iteration to be run
loop = 2;

% number of poses to be simulated
% poses = 120;
poses = 24;

%% Main code

% generate ground truth data
[subA, supA, subB, supB, X, Y] = funPreparingAXXB(poses);

for m = 1:3
    m
    % phantom model
    [model] = funPhantomModel(m);

    
    for n = 1:loop
        % out of plane error (mid-plane detection is required in pose 1)
        errorOOP = MagErrorOOP*randn(1,1);
        
        p = zeros(4,npoint, poses);
        p2d = zeros(3,npoint, poses);
        nsp = zeros(npoint, poses);
        for i = 1:poses
            for point = 1:npoint
                
                % defining ground truth position of phantom from ultrasound
                % image frame
                p(:,point,i) = subA(:,:,i)*model(:,point);
                if p(2,point,i) < 0
                    p(2,point,i)
                    pause(100)
                end
                % point position appeared on 2D ultrasound image view
                p2d(1,point,i) = p(1,point,i);
                p2d(2,point,i) = sqrt(p(2,point,i)^2+p(3,point,i)^2);
                p2d(3,point,i) = 0;
                nsp(point,i) = p(3,point,i);
            end
        end
        
        % segmented point position
        sp(1:2,:,1:eposes) = p2d(1:2,:,1:eposes)+MagErrorSeg*rand(2,npoint,eposes);
        sp(3,:,1:eposes) = p2d(3,:,1:eposes);

        % optimization setting
        opts = optimset('display','off','MaxFunEvals',5000,'MaxIter',1000);
        
        % initialization
        residual = 1;
        initV = [zeros(6,1);nsp(1:npoint,1)];
        [initV(3), initV(2), initV(1)] = decompose_rotation_d(subA(:,:,1));
        initV(4:6) = subA(1:3,4,1);
        minCost = inf;
        
        % solving full A for the first pose
        for l = 1:50
            [initV,residual] = lsqnonlin(@funRecoverPartialA_OneMidPlane,initV,[],[],opts,model(:,1:npoint),sp(:,:,1),errorOOP);
            initV(6+1:6+npoint) = [nsp(1:npoint,1) + real(10*randn(npoint,1))];
            if residual < minCost
                minCost = residual;
                v = initV;
            end
        end
        
        fullA1 = (buildT(v(1:6)));
        partA = zeros(4,4,eposes);
        rot = zeros(4,4,eposes);
        x = zeros(1,eposes);
        y = zeros(1,eposes);
        z = zeros(1,eposes);
        
        for i = 2:eposes
            % initialization
            residual = 1;
            initV = [zeros(6,1);nsp(1:npoint,i)];
            [initV(3), initV(2), initV(1)] = decompose_rotation_d(subA(:,:,i));
            initV(4:6) = subA(1:3,4,i);
            minCost = inf;
            count = 0;
            
            % solving for partial As for the rest of poses
            for l = 1:5
                count = count + 1;
                [tempV,residual] = lsqnonlin(@funRecoverPartialA,initV,[],[],opts,model(:,1:npoint),sp(:,:,i));
                initV(6+1:6+npoint) = [nsp(1:npoint,i) + real(count*randn(npoint,1))];
                residual;
                if residual < minCost
                    minCost = residual;
                    v = tempV;
                end
            end
            partA(:,:,i) = (buildT(v(1:6)));
            rot(:,:,i) = (partA(:,:,i))*inv(subA(:,:,i));
            [x(i), y(i), z(i)] = decompose_rotation_d(rot(:,:,i));
        end
        
        % recovering the unknown rotations of As using one full A and
        % partically recovered As
        
        % initialization
        opts=  optimset('display','on','MaxFunEvals',50000,'MaxIter',10000);
        initV2 = [0 x(2:eposes)]';
        
        % set the first pose to be fully recovered A
        partA(:,:,1) = fullA1;
        minCost = inf;
        
        % recover the unknown rotation in A
        count = 0;
        for l = 1
            count = count + 1;
            [tempV,cost,residual] = lsqnonlin(@funRecoverA_2invariants,initV2,[],[],opts,partA(:,:,1:eposes),subB(:,:,1:eposes),fullA1);
            if residual < minCost
                initV2 = [0 x(2:eposes)]' + real(round(count/3)*randn(eposes,1));
                minCost = cost;
                v2 = tempV;
            end
        end
        
%         toc
        RotA = zeros(4,4,eposes);
        newA = zeros(4,4,eposes);
        
        % recover full transformation of As
        for i = 1:eposes
            RotA(:,:,i) = buildT([0 0 v2(i) 0 0 0]);
            newA(:,:,i) = inv((partA(:,:,i))\RotA(:,:,i));
        end
        newA(:,:,1) = fullA1(:,:);
        
        % reconstruct superscript As using recovered subscript As
        count = 0;
        for i = 1:eposes
            for j = 1:eposes
                if i ~= j
                    count = count + 1;
                    rsupA(:,:,count) = (newA(:,:,i))*inv(newA(:,:,j));
                    esupB(:,:,count) = inv(subB(:,:,i))*subB(:,:,j);
                end
            end
        end
        
        % solving X through AXXB
        [ reX, residual ] = solveAXXB(rsupA(:,:,1:count), esupB(:,:,1:count));

        % setting the number of poses to be tested
        tposes = poses-eposes;
        
        % setting ground truth point
        c = subB(:,:,poses-tposes+1)*inv(X)*p(:,1,poses-tposes+1);
        
        % initialization reconstruction precision and accuracy
        rp = zeros(4,tposes);
        accup = zeros(4,tposes);
        
        % points that used to compute reconstruction precision (rp) and accuracy (accup)
        for i = 1:tposes
            rp(:,i) = subB(:,:,poses-tposes+i)*inv(reX)*p(:,1,poses-tposes+i);
            accup(:,i) = abs(rp(:,i)-c);
        end

        % computing accuracy
        t_accu_mean = mean(accup(1:3,:),2);
        accu_mean(:,n,m) = norm(t_accu_mean);
        
        % computing reconstruction precision
        x_std = std(rp(1,:));
        y_std = std(rp(2,:));
        z_std = std(rp(3,:));
        rp_std(:,n,m) = norm([x_std y_std z_std]);
    end
    
    % averaged reconstruction precision and accuracy with multiple trials
    precision(:,m) = mean(rp_std(:,:,m))
    accuracy(:,m) = mean(accu_mean(:,:,m))
end