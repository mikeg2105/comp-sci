%To review the results load the matlab file

for ic=1:1:4
    subplot(2,2,ic);
    surf(res{ic}, 'LineStyle', 'none');
%    xlabel('shift');
%    ylabel('tscale');
end