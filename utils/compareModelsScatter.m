function compareModelsScatter(folder1, folder2, title1, title2)

% create the filenames for each folder
filename1 = fullfile(folder1, 'correlation-Val.nii');
filename2 = fullfile(folder2, 'correlation-Val.nii');

% load the data for both files
hdr1 = spm_vol(filename1);
hdr2 = spm_vol(filename2);
data1 = spm_read_vols(hdr1);
data2 = spm_read_vols(hdr2);

% plot the scatter plot
figure;
scatter(data1(:), data2(:));
xlabel(title1);
ylabel(title2);
xlimits = xlim();
ylimits = ylim();
newlimits = [min(xlimits(1),ylimits(1)) max(xlimits(2),ylimits(2))];
xlim(newlimits);
ylim(newlimits);
line;