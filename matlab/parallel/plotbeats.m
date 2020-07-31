 
for ic=1:4
    infile= [ 'beat' ,num2str(ic) , '.mat' ];
    load(infile);
    subplot(2,2,ic);
    surf(wave, 'LineStyle', 'none');
    xlabel([ 'shift' , num2str(ic) ]);
    ylabel('tscale');
    
end
