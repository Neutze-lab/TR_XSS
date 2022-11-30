% MatWAXS interpolator
clear all
close all

addpath('/path/to/inputs');
addpath('/path');
%% SMX open/closed SFX ground + cis retinal - aligned in space and all have trans retinal identical x y z 
number_of_iterations = 1;
output_prefix =['name_of_intermediate'];
exagerate = [number1,number2];
mkdir(output_prefix)
cd(output_prefix)
resting_state = '5b6z_resting_aligned.pdb';
exci_state = '5b6z_excited_aligned.pdb';
new_ground = '5b6z_resting_aligned.pdb';
output_dir = pwd;

%% Residue range to exagerate
res_range = [77;91;152;186] %new range 02/2022
%res_range = [81;91;151;181]; %Helix C & Helix E-F
% res_range = [81;91]; %Helix C
% res_range = [143;224] %Helix EFG



% residues = [81,91; 151,181]
% save('residues.mat','residues');
x = zeros(size(exagerate,2),1);

    for n = 0:number_of_iterations       
        
    step_shift = exagerate/number_of_iterations;
    x = step_shift*n;
        if (n == 0)
            x = zeros(size(exagerate,2),1);
        end
    fprintf('%d %d %d\n',n,x)
    output_file = sprintf('%s_%d.pdb',output_prefix,n);
    logfile = output_file(1:end-4);

       vectoradd_v3(x,char(resting_state), char(exci_state), output_file(1:end-4), output_dir, res_range, new_ground);   
        
    end

% remove extra file and rename     
delete *0.pdb
movefile(output_file,strcat(output_prefix,'.pdb'))
cd ..

function new_pos_av = vectoradd_v3(m,avst,extst,output,output_dir,res_range,new_ground)

if (~exist(new_ground,'file') && string(new_ground) ~= "none")
    sprintf('Your different ground file does not exist, please write "none" if you wish to not use the different ground specfication');
    return
end
    
n = 0;
pos_av = [];
pos_st1 = [];
vec1 = [];

mod_av = pdbreader(avst);
mod_av2 = pdbreader('6rqp_SMX_closed_aligned.pdb');
mod_av = mod_av(ismember(mod_av(:,1),"ATOM"),:);
mod_av2 = mod_av2(ismember(mod_av2(:,1),"ATOM"),:);
mod_st1 = pdbreader(extst);
mod_st1_2 = pdbreader('6rph_SMX_open_aligned_NOH_LOOP.pdb');
mod_st1 = mod_st1(ismember(mod_st1(:,1),"ATOM"),:);
mod_st1_2 = mod_st1_2(ismember(mod_st1_2(:,1),"ATOM"),:);

mod_temp = mod_av;

% [mod_av,mod_st1] = pdb_unfucker(mod_av,mod_st1);

if (exist(new_ground,'file'))
    diff_ground_mod_av = pdbreader(new_ground);
    diff_ground_mod_av = diff_ground_mod_av(ismember(diff_ground_mod_av(:,1),"ATOM"),:);
%     [mod_av,diff_ground_mod_av] = pdb_unfucker(mod_av,diff_ground_mod_av);
%     [mod_st1,diff_ground_mod_av] = pdb_unfucker(mod_st1,diff_ground_mod_av);
end

Ca_resting = find(mod_av(:,3)=='CA');
Ca_ground = find(mod_st1(:,3)=='CA');

pos_av = double(mod_av(:,get_cartesian(mod_av)));
pos_av2 = double(mod_av2(:,get_cartesian(mod_av2)));
pos_st1 = double(mod_st1(:,get_cartesian(mod_st1)));
pos_st1_2 = double(mod_st1_2(:,get_cartesian(mod_st1_2)));
res_group = double(mod_av(:,5));
    

% BIT WHERE WE CREATE NEW VECTOR WITH EXAGERATIONS

curr_res_range1 = res_range(1);
curr_res_range2 = res_range(2);
exag_iterator = 1;
exag = m(exag_iterator);
res_point = 2;

    
    for i = 1:length(res_group)
            if (res_group(i) <= curr_res_range1 || res_group(i) >= curr_res_range2)
                vec1(i,1:3) = (pos_av(i,:)-pos_av(i,:)).*exag;
            else
                vec1(i,1:3) = (pos_st1(i,:)-pos_av(i,:)).*exag;
            end            

            if (size(res_range,1)>res_point)
                if (res_group(i) >= curr_res_range2)
                    curr_res_range1 = res_range(res_point+1);
                    curr_res_range2 = res_range(res_point+2);
                    pos_av = pos_av2;
                    pos_st1 = pos_st1_2;
                    exag_iterator = exag_iterator + 1;
                    exag = m(exag_iterator);
                    res_point = 2 + 2;
                end 
            end
    end


if (exist(new_ground,'file'))
    output_pdb = pdbread(new_ground);
    diff_ground_pos_av = double(diff_ground_mod_av(:,get_cartesian(diff_ground_mod_av)));
    new_pos_av = round(diff_ground_pos_av + vec1,3);
    diff_ground_mod_av(:,6:8) = new_pos_av;
    new_mod_av = diff_ground_mod_av;
else
    output_pdb = pdbread(avst);
    new_pos_av = round(pos_av + vec1,3);
    mod_av(:,6:8) = new_pos_av;
    new_mod_av = mod_av;
end

% BIT WHERE WE SAVE THE NEW VECTOR

for a = 1:size(new_pos_av,1)
    output_pdb.Model.Atom(a).X = double(new_mod_av(a,6));
    output_pdb.Model.Atom(a).Y = double(new_mod_av(a,7));
    output_pdb.Model.Atom(a).Z = double(new_mod_av(a,8));

end

output = cat(2,output,'.pdb');
pdbwrite(output,output_pdb);



function [ground_atoms,excited_atoms] = pdb_unfucker(ground_atoms, excited_atoms)
    residue_number_ground = unique(ground_atoms(:,5),'rows');
    residue_number_excited = unique(excited_atoms(:,5),'rows');    
    ground_atoms = remove_residues(residue_number_ground,residue_number_excited,ground_atoms);
    excited_atoms = remove_residues(residue_number_excited,residue_number_ground,excited_atoms);
    residue_number_ground = unique(ground_atoms(:,5),'rows');
    residue_number_excited = unique(excited_atoms(:,5),'rows');       
    tic
    
    for i = 1:size(residue_number_ground,1)
        current_res_num = residue_number_ground(i);
        res_atoms_ground = ground_atoms(:,3);
        res_atoms_excited = excited_atoms(:,3);     
            ground_index = find(ground_atoms(:,5) == residue_number_ground(i));
            excited_index = find(excited_atoms(:,5) == residue_number_excited(i));            
            residue_elements_ground = res_atoms_ground(ground_index);
            residue_elements_excited = res_atoms_excited(excited_index);         
            ground_atoms = remove_atoms(residue_elements_ground,residue_elements_excited,ground_atoms,ground_index);
            excited_atoms = remove_atoms(residue_elements_excited,residue_elements_ground,excited_atoms,excited_index);
    end
    [ground_atoms,excited_atoms] = order_atoms(ground_atoms,excited_atoms,current_res_num);
    toc
    
end    
     
function atom_set = remove_residues(ground_residues,excited_residues,atom_set)
    diff_elements = setdiff(ground_residues,excited_residues);
        for a = 1:size(diff_elements,1) 
            index_remove = find(diff_elements(a) == atom_set(:,5));
            atom_set(index_remove,:) = [];
        end 
end

function atom_set = remove_atoms(remove_these_atoms,keep_these_atoms,atom_set,atom_index)
    atoms_to_remove = setdiff(remove_these_atoms,keep_these_atoms);
    if(size(atoms_to_remove,1) ~= 0) 
%         for b = 1:size(atoms_to_remove,1);
            [~,index_remove] = intersect(remove_these_atoms,atoms_to_remove);
            atom_set(atom_index(index_remove),:) = [];
%         end
    end
end

function [atom_set1,atom_set2] = order_atoms(atom_set1,atom_set2,res_num)
    res_nums1 = atom_set1(:,5);
    res_nums2 = atom_set2(:,5);
    index_res_num1 = find(res_num == res_nums1);
    index_res_num2 = find(res_num == res_nums2);
    atom_identifiers1 = atom_set1(index_res_num1,:);
    atom_identifiers2 = atom_set2(index_res_num2,:);
    atom_identifiers1 = sortrows(atom_identifiers1,3);
    atom_identifiers2 = sortrows(atom_identifiers2,3);
    atom_set1(index_res_num1,:) = atom_identifiers1;
    atom_set2(index_res_num2,:) = atom_identifiers2;
end
    
end


function output = get_cartesian(data)

% file = '/path';
% file2 = '/path';
% data = get_mod_dat(file);

a = double(data(1,:));
a(isnan(a)) = 1;
b = size(a,2);
c = [];

 for d = 1:b
    if(mod(a(d),1)~=0)
        c = [c;d];
        if(size(c,1)==3)
            output = c;
        end
    end
 end
end

function rhysed_mat = pdbreader(input_path)
% fprintf('\n%s',input_path);
input = string(fileread(input_path));
first_pass(:,1) = regexp(input,"(?m)^(ATOM  |CONECT|HETATM|TER.{0,3}|END\s{0,3})(.{5})?(.{5})?(.{4})?(.{2})?(.{4})?(.{12})?(.{8})?(.{8})?(.{6})?(.{6})?(.{12})?.*",'tokens','dotexceptnewline');
rhysed_mat = strtrim(vertcat(first_pass{:,1}));
rhysed_mat(:,5) = []; 
end



