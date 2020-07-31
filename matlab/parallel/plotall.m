function  plotall ( ii ) 
% expects an integer argument. Reads wave pattern from a series of files
% named  wave(number).mat  and plots them one at a time.
%
if nargin ~= 1
    error( 'plotall:IncorrectNumberInputs', ' An integer number expected' )
    return
end
for i=1:ii
    ichar = num2str( i ) ;
    infile= ['wave' ichar '.mat' ];
    load ( infile );
    surf ( wave1 , 'LineStyle', 'none'); 
    title ( [ 'Plotting series: ' , ichar ] )
    shg
    pause(2)
end
 