function void = process_intensity_files(cntmax,prefix,suffix,cut)
%   function a_sorted_cell = sort_cell_file_names(a_cell)
%   Holiday fun would be to come back and write this script with structures.

fprintf('Process intensities has started for %s\n',suffix)
all_intensities = [];
for j = 1:cntmax 
    if (strcmp(cut,"nope"))
    fprintf('nope has worked\n')
    search_reg_expr=prefix+"*_"+string(j)+"*_"+suffix+"*[0-9].int";
    else
    search_reg_expr=prefix+"*_"+string(j)+"*_"+suffix+"*_noplace*SphereCut.int";
    end
    [~,~] = unix(sprintf('find . -maxdepth 1 -name "%s" > matlab_filenames.list',search_reg_expr));
    fid = fopen('matlab_filenames.list')
    file_names = textscan(fid,'%s');
    fclose(fid);
%    delete('matlab_filenames.list')
    files = file_names{1,1};
    files = sort(string(files));

    dim1 = size(dlmread(files{1,1}," ",1,0),1);
    run_intensities = zeros(dim1,4,size(files,1));
    for m = 1:size(files,1)
        fprintf('%s\n',files{m,1})
        a = dlmread(files{m,1}," ",1,0);
        run_intensities(:,:,m) = a(:,[3 5 7 9]);
    end
    all_intensities = cat(3,all_intensities,run_intensities);
end


q = all_intensities(:,1,1);
intSol = squeeze(all_intensities(:,2,:));
intVac = squeeze(all_intensities(:,3,:));
solvScatt = squeeze(all_intensities(:,4,:));

if (strcmp(cut,"nope"))
save(char(prefix+"_intensities_"+suffix+".mat"),'q','intSol','intVac','solvScatt')
else
save(char(prefix+"_intensities_"+suffix+"_noplaceÅ_SphereCut.mat"),'q','intSol','intVac')
end

    
end    

