clc; close all; clear;
%%


S = dir(fullfile('results/','*.mat')); 
%
for k = 1:numel(S)
    N = convertCharsToStrings(S(k).name);
    if contains(N,"struct")
        load(N);
    else
        continue;
    end
    
    N_split = split(N,'_');
    scene = N_split(end-1);
    
    method = split(N_split(end),'.');
    method = method(1);
    
    field_names_cell = fieldnames(meta_results_struct);
    for ii = 1:length(field_names_cell)
        s_id = field_names_cell{ii};
        
        temp = getfield(meta_results_struct,s_id);
        
%         disp(join([s_id " "...
%             num2str(temp.start_bits/temp.end_bits_scores)]));
        
        
    end
    
    plot_meta_results_func_v2(meta_results_struct, scene, method)
    
end