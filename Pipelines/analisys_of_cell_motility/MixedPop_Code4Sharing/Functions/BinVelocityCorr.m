function [V_binned, V_binned_err, V_dR] = BinVelocityCorr(n_frames,R_bins,delta_R,cell_Rnm,cell_VnVm)

V_binned = NaN(numel(R_bins),n_frames-1);
V_binned_err = NaN(numel(R_bins),n_frames-1);

V_dR = cell(numel(R_bins),1);

for jj = 1 : n_frames-1

    for hh = 1 : numel(R_bins)

        idx_R = cell_Rnm{jj} < (R_bins(hh) + delta_R) & cell_Rnm{jj} > (R_bins(hh) - delta_R);

        V_dR{hh} = cell_VnVm{jj}(idx_R);

        V_binned(hh,jj) = mean(cell_VnVm{jj}(idx_R),'omitnan');

        V_binned_err(hh,jj) = std(cell_VnVm{jj}(idx_R),0,'all','omitnan');

    end

end

end