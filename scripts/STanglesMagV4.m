
function [angleMagTuned] = STanglesMagV4(spatial_tunedROI_A,spatial_selectivity_A,spatial_tunedROI_B,spatial_selectivity_B,singleCtxFlag)

%identify which cells retained spatial tuning, gained, lost and never got tuned
%calculate the magnitude of spatial vectors before and after
%calculate the tuning angle relative to 0 deg being +x axis
%calculate the centroid shift of cells that retained tuning
%comparisons relative to the first 'context'/set of laps
%context A or B can be integer 1, 2, 3
%calculate the magnitude and angles

%spatial_tunedROI = binary vector mask
%spatial_selectivity = complex double vectors

if singleCtxFlag == 1
    %get idxs of the spatially tuned cells
    tunedROIidx = find(spatial_tunedROI_A == 1);
    
    for i=1:length(tunedROIidx)
        rX = real(spatial_selectivity_A(i));
        rY = imag(spatial_selectivity_A(i));
        tuneAnglePrior(i) = atan2d(rY,rX);
        if (tuneAnglePrior(i) < 0)
            tuneAnglePrior(i) = tuneAnglePrior(i) + 360;
        end
        
      magPre(i) = abs(spatial_selectivity_A(i));  
    end
    
   angleMagTuned.tuneAngle = tuneAnglePrior;
   angleMagTuned.mag = magPre;

elseif (singleCtxFlag == 0)
    
    retainTuning = find(spatial_tunedROI_A + spatial_tunedROI_B == 2);
    newTuned = find(spatial_tunedROI_B*2 - spatial_tunedROI_A == 2);
    loseTuned = find(spatial_tunedROI_B - spatial_tunedROI_A == -1);
    neverTuned = find(spatial_tunedROI_A + spatial_tunedROI_B == 0); 
    
    %get mask of tuned cells
    tunedROIs = {retainTuning, newTuned, loseTuned, neverTuned};

    %calculate the angle before and after exposure to odor for cells that
    %retain spatial tuning
    %preAngles = {};
    
    for ii = 1:length(tunedROIs)
        tuneAnglePrior = [];
        tuneAnglePost = [];
        
        for i=1:length(tunedROIs{ii})
            if (ii == 1 || ii == 3)
                rX = real(spatial_selectivity_A(tunedROIs{ii}(i)));
                rY = imag(spatial_selectivity_A(tunedROIs{ii}(i)));
                tuneAnglePrior(i) = atan2d(rY,rX);
                if (tuneAnglePrior(i) < 0)
                    tuneAnglePrior(i) = tuneAnglePrior(i) + 360;
                end
            end
            %figure; compass(spatial_selectivity{1,1}(retainTuning(i)))
            
            if (ii == 1 || ii == 2)
                rX = real(spatial_selectivity_B(tunedROIs{ii}(i)));
                rY = imag(spatial_selectivity_B(tunedROIs{ii}(i)));
                tuneAnglePost(i) = atan2d(rY,rX);
                if (tuneAnglePost(i) < 0)
                    tuneAnglePost(i) = tuneAnglePost(i) + 360;
                end
            end
            %figure; compass(spatial_selectivity{1,3}(retainTuning(i)))
            
        end
        
        if (ii == 1)
            angleDiffPP = tuneAnglePost - tuneAnglePrior;
        end
        
        if (ii == 1 || ii == 3)
            preAngles{ii} = tuneAnglePrior;
        end
        
        if (ii == 1 || ii == 2)
            postAngles{ii} = tuneAnglePost;
        end
        
    end
    
    %calculate the magnitudes
    for ii = 1:length(tunedROIs)
        for i=1:length(tunedROIs{ii})
            if (ii == 1 || ii == 3) %tuned and pre-tuned only
                magPre{ii}(i) = abs(spatial_selectivity_A(tunedROIs{ii}(i)));
            end
            
            if (ii == 1 || ii == 2) %tuned and post-tuned only
                magPost{ii}(i) = abs(spatial_selectivity_B(tunedROIs{ii}(i)));
            end
            
        end
    end
    
    %package variables into struct
    angleMagTuned.tunedROIs = tunedROIs;
    angleMagTuned.angleDiffPP = angleDiffPP;
    angleMagTuned.preAngles = preAngles;
    angleMagTuned.postAngles = postAngles;
    angleMagTuned.magPre = magPre;
    angleMagTuned.magPost = magPost;
    
end
    
end