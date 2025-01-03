function [x1,x2, desc_X, fA, dA, fB, dB] = feature_extraction(im1,im2)
        
        % Return a homogenous projections

        % Extract the features using SIFT
        [fA dA] = vl_sift(single(rgb2gray(im1)) );
        [fB dB] = vl_sift(single(rgb2gray(im2)) );
        
        % Compute the matches
        matches = vl_ubcmatch(dA, dB);
        
        % Get the 2D points correspondences on each images' pair
        xA = fA(1:2 , matches(1 ,:));
        xB = fB(1:2 , matches(2 ,:));
        
        % Convert the vectors to be homogenous for this pair of images 
        % Add a homogenous dimension to the points
        x1 = [xA; ones(1,length(xA))];
        x2 = [xB; ones(1,length(xB))];
    
        % Save descriptors for matched points (3D points will correspond to these)
        % Extract descriptors only for matched points in image 1
        desc_X = dA(:, matches(1, :)); 
       
end