function voxIdxs = wrapperVoxelSelection(options, subject)

maskFolder = fullfile(options.maskFolder, subject);
if strcmp(options.voxelSelectionMask,'')
    % create the folder name where the brain mask is expected to be
    maskWildcard = fullfile(maskFolder, 'voxelSelectionMask.nii');
    maskFilenames = dir(maskWildcard);
    assert(numel(maskFilenames)==1);
    maskFilename = fullfile(maskFolder, maskFilenames(1).name);
else
    maskFilename = fullfile(maskFolder, options.voxelSelectionMask);
end

% load the mask data
header = spm_vol(maskFilename);
maskData = spm_read_vols(header);

% now find the indices that have a 1 indicating the brain exists in that
% voxel
voxIdxs = find(maskData==1);