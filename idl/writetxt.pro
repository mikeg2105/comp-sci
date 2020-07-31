

args=command_line_args()
index=long(args(0))
print,'index is ',index
print,'args is ',args

filename=args+'.out'

openw,1,filename

printf,1,'This is the index',index
printf,1,'This is the string',args

close,1

end

