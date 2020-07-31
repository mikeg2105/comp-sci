

args=command_line_args()
index=long(args(0))
print,'index is ',args

;Read data from SDSS spectra
;wget "http://www.voservices.net/spectrum/search_details.aspx?format=csv&id=ivo%3a%2f%2fjhu%2fsdss%2fdr6%2fspec%2f2.5%23XXSPECOBJIDxx"

;************************************************
;
;specsedscript file contents
;s/,/ /g
;************************************************


;specid='96486455647404032'
specid='587104588549586944'

ssidar=strarr(7)
ssidar(0)='96486455647404032'
ssidar(1)='102398707621691392'
ssidar(2)='79034666178314240'
ssidar(3)='79880582720389120'
ssidar(4)='80443407029239808'
ssidar(5)='96486453780938752'
ssidar(6)='84102689763360768'

for i=5,6 do begin
;i=3
specid=ssidar(i)
wgetcommand='wget "http://www.voservices.net/spectrum/search_details.aspx?format=csv&id=ivo%3a%2f%2fjhu%2fsdss%2fdr6%2fspec%2f2.5%23'+specid+'"'
mvcommand='mv search_details*'+specid+' temp.csv'
cpcommand2='cp '+specid+'.csv temp.csv'
cpcommand='cp testspec.csv '+specid+'.csv'
sedcommand='sed -f specsedscript < temp.csv > testspec.csv'
spawn,wgetcommand
spawn,mvcommand
;spawn,cpcommand2
spawn,sedcommand
;spawn,'mv temp.csv testspec.csv'
spawn,cpcommand

endfor

end

