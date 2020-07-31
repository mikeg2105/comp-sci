%To review the results load the matlab file

for ic=1:3:11
    subplot(2,2,ic);
    surf(res{ic}, 'LineStyle', 'none');
%    xlabel('shift');
%    ylabel('tscale');
end