function [x1,x2, desc_X, fA, dA, fB, dB, matches] = feature_extraction(im1,im2)
    % FEATURE_EXTRACTION - Extracts keypoints, descriptors, and matches between two images.
    %
    % Inputs:
    %   im1, im2 - Input images (grayscale or RGB).
    %
    % Outputs:
    %   x1, x2   - 2D coordinates of matched keypoints in both images (2xN).
    %   desc_X   - Descriptors of matched keypoints in the first image (128xN).
    %   fA, dA   - Keypoints and descriptors for the first image.
    %   fB, dB   - Keypoints and descriptors for the second image.

        
        % Return a 2D corresponces and other ....

        % Extract the features using SIFT
        % 1 specifies a threshold for rejecting low-contrast keypoints detected by SIFT.

        [fA, dA] = vl_sift(single(rgb2gray(im1)),'PeakThresh', 2);
        [fB, dB] = vl_sift(single(rgb2gray(im2)),'PeakThresh', 2);
        
        % Compute the matches
        matches = vl_ubcmatch(dA, dB);
        
        % Get the 2D points correspondences on each images' pair
        x1 = fA(1:2 , matches(1 ,:));
        x2 = fB(1:2 , matches(2 ,:));
  
    
        % Save descriptors for matched points (3D points will correspond to these)
        % Extract descriptors only for matched points in image 1
        desc_X = dA(:, matches(1, :)); 
       
end