function [] = plot_meta_results_func_v2(meta_results_struct, scene, dr_method)

% Save meta results struct
warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir results;


fig_handle = figure;
hold on;

field_names_cell = fieldnames(meta_results_struct);
for ii = 1:length(field_names_cell)
    s_id = field_names_cell{ii};

    % unpack struct
    qi_load = meta_results_struct.(s_id).("qi_load");
    qi_comp = meta_results_struct.(s_id).("qi_comp");
    qi_load_res = meta_results_struct.(s_id).("qi_load_res");
    qi_comp_res = meta_results_struct.(s_id).("qi_comp_res");
    
    % scatter plot
    sz = 20;
    mkr = 'o';

    field_name = meta_results_struct.(s_id).wname + " " + ...
        num2str(meta_results_struct.(s_id).wtreshold);
    
    clr_load = 'g';
    clr_comp = 'r';
    clr_load_res = 'c';
    clr_comp_res = 'b';
    
    cr_load = 100/meta_results_struct.(s_id).("dr_loadings"); % Change
    
    cr_comp = meta_results_struct.(s_id).start_bits / ...
                  (meta_results_struct.(s_id).end_bits_scores);
    
	cr_load_res = meta_results_struct.(s_id).start_bits / ...
                  (meta_results_struct.(s_id).end_bits_res + ...
                  meta_results_struct.(s_id).end_bits_scores);
    
    cr_comp_res = 2*100/meta_results_struct.(s_id).("dr_loadings"); % Change
    
    
    s_load = scatter(qi_load.psnr,cr_load,...
        sz,clr_load,mkr);
    s_comp = scatter(qi_comp.psnr,cr_comp,...
        sz,clr_comp,mkr);
    s_load_res = scatter(qi_load_res.psnr, cr_load_res, ...
        sz,clr_load_res,mkr);
    s_comp_res = scatter(qi_comp_res.psnr,cr_comp_res,...
        sz,clr_comp_res,mkr);
    
    row = dataTipTextRow(field_name,0);
    
    s_load.DataTipTemplate.DataTipRows(1) = row;
    s_load.DataTipTemplate.DataTipRows(2).Label = 'CR';

    s_comp.DataTipTemplate.DataTipRows(1) = row;
    s_comp.DataTipTemplate.DataTipRows(2).Label = 'CR';

    s_load_res.DataTipTemplate.DataTipRows(1) = row;
    s_load_res.DataTipTemplate.DataTipRows(2).Label = 'CR';

    s_comp_res.DataTipTemplate.DataTipRows(1) = row;
    s_comp_res.DataTipTemplate.DataTipRows(2).Label = 'CR';
    
    fields_qi = fieldnames(qi_load);
    for jj = 1:length(fields_qi)
        field_name = fields_qi{jj};
        
        row = dataTipTextRow(field_name,qi_load.(field_name));
        s_load.DataTipTemplate.DataTipRows(end+1) = row;
        
        row = dataTipTextRow(field_name,qi_comp.(field_name));
        s_comp.DataTipTemplate.DataTipRows(end+1) = row;
        
        row = dataTipTextRow(field_name,qi_load_res.(field_name));
        s_load_res.DataTipTemplate.DataTipRows(end+1) = row;
        
        row_res = dataTipTextRow(field_name,qi_comp_res.(field_name));
        s_comp_res.DataTipTemplate.DataTipRows(end+1) = row_res;
    end
    
end    

title(join([scene " - " dr_method " - True Lossy Compression of Scores"]));
xlabel("Mean Bandwidth Peak Signal to Noise Ratio");
ylabel("Compression Ratio");

set(fig_handle,'position',[200, 300, 1300, 1000])
grid on;
hold off;


fname_figure = "results/meta_results_"+date+"_"+scene+"_"+dr_method+".fig";
saveas(fig_handle, fname_figure);

fname_final = "figs/meta_results_"+scene+"_"+dr_method+".pdf";
exportgraphics(fig_handle,fname_final)


end
