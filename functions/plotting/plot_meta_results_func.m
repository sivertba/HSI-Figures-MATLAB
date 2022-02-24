function [] = plot_meta_results_func(meta_results_struct, scene, dr_method)

% Save meta results struct
warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir results;

fname_figure= "results/meta_results_struct_"+scene+"_"+dr_method+".mat";

save(fname_figure,...
    'meta_results_struct')

fig_handle = figure;
hold on;

for ii = 1:numel(fieldnames(meta_results_struct))
    disp(ii);
    % unpack struct
    qi = meta_results_struct.("id_"+num2str(ii)).("qi");
    qi_res = meta_results_struct.("id_"+num2str(ii)).("qi_res");
%     ar = meta_results_struct.("id_"+num2str(ii)).("ar");
    cp = meta_results_struct.("id_"+num2str(ii)).("cp");
    
    cr = 1/meta_results_struct.("id_"+num2str(ii)).ar.compressedSizewoRes_comp;
    psnr = qi.psnr;

    cr_res = 1/meta_results_struct.("id_"+num2str(ii)).ar.compressedSize_comp;
    psnr_res = qi_res.psnr;
    
    
    % scatter plot
    sz = 20;
    mkr = 'o';
    fields_cp = fieldnames(cp);

       
    if strcmp(cp.method, "wavelet")
        clr = 'r';
        field_name = cp.method + " " + cp.wname;
    elseif strcmp(cp.method, "fft")
        clr = 'g';
        field_name = cp.method;
    elseif strcmp(cp.method, "None")
        cr = 1/meta_results_struct.("id_"+num2str(ii)).ar.compressedSizewoRes;
        cr_res = 1/meta_results_struct.("id_"+num2str(ii)).ar.compressedSize;
        clr = 'b';
        field_name = cp.method;
    else
        warning("Uknown method attempted to be plotted");
        clr = 'k';
    end
    
    s = scatter(psnr,cr,sz,clr,mkr);
    s_res = scatter(psnr_res,cr_res,sz,clr,mkr,'filled');
    
    row = dataTipTextRow(field_name,0);
    s.DataTipTemplate.DataTipRows(1) = row;
    s.DataTipTemplate.DataTipRows(2).Label = 'CR';

    s_res.DataTipTemplate.DataTipRows(1) = row;
    s_res.DataTipTemplate.DataTipRows(2).Label = 'CR';

    for jj = 1:length(fields_cp)
        field_name = fields_cp{jj};        
        val = cp.(field_name); 
        if ischar(val)
            continue
        end
        row = dataTipTextRow(field_name,val);
        disp(row);
        s.DataTipTemplate.DataTipRows(end+1) = row;
        s_res.DataTipTemplate.DataTipRows(end+1) = row;
    end
    
    fields_qi = fieldnames(qi);
    for jj = 1:length(fields_qi)
        field_name = fields_qi{jj};    
        val = qi.(field_name); 
        val_res = qi.(field_name); 
        if ischar(val)
            continue
        end
        row = dataTipTextRow(field_name,val);
        row_res = dataTipTextRow(field_name,val_res);
        s.DataTipTemplate.DataTipRows(end+1) = row;
        s_res.DataTipTemplate.DataTipRows(end+1) = row_res;
    end
    row = dataTipTextRow("residuals",0);
    s.DataTipTemplate.DataTipRows(end+1) = row;
    row_res = dataTipTextRow("residuals",1);
    s_res.DataTipTemplate.DataTipRows(end+1) = row_res;
end    

title(join([scene " - " dr_method " - Lossy Compression of Scores"]));
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
