%  recovering variables from a cell array .... 
for ic=1:4
     
    subplot(2,2,ic);
    surf(wave{ic} , 'LineStyle', 'none');
    xlabel([ 'shift' , num2str(ic) ]);
    ylabel('tscale');
    
end
