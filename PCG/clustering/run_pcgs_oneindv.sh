#!/bin/bash

indv=$1
binnum=$2
target_gene=$3

echo $target_gene

rm -rf   ${indv}"_bin"${binnum}
mkdir -p ${indv}"_bin"${binnum}
cd       ${indv}"_bin"${binnum}
#pwd


#target_gene='gene_name_xu230330_gen_data.txt'
pos_file="cell_pos_db_filter_singlecell/${indv}.txt"

awk -v binnum=$binnum 'BEGIN{
        yxmin=100;
        yxmax=0;
     }
     {
        if(y<10 &&y>-10&&$1<10)
           yxmin=$1;  
        if($1>max) 
           max=$1;  
        if(y<10 &&y>-10&&$1>yxmax)yxmax=$1; 
     }
     END{
        printf("xmin=%d\n",yxmin); 
        printf("xmax=%d\n",yxmax);
        printf("binsize=%s\n",(yxmax-yxmin+1)/binnum);
     }' $pos_file  > ${indv}.bininfo.txt

xmin=`grep xmin ${indv}.bininfo.txt | awk -F "=" '{print $2}'`
xmax=`grep xmax ${indv}.bininfo.txt | awk -F "=" '{print $2}'`
binsize=`grep binsize ${indv}.bininfo.txt | awk -F "=" '{print $2}'`

awk -v xmin=$xin -v xmax=$max -v binsize=$binsize -v binnum=$binnum   '{
        binid = ($1-xmin)/binsize;
        if(binid>=0||binid<binnum){
            t[int(binid)]+=1;
        }
    }
    END{
        for(x=0;x<binnum;x++)
            if(x in t)
                print x,t[x];
            else
                print x,0;
    }' $pos_file > ${indv}.cells_in_bin.txt


while read -r smes
do
    smes_file="${indv}/${smes}.txt"
    #echo $smes_file
    if [[ -e $smes_file ]] ; then
        awk -v xmin=$xin -v xmax=$max -v binsize=$binsize -v gene=$smes -v binnum=$binnum '
            BEGIN{
               print gene;
            }
            {
                if(FILENAME==ARGV[1]) {
                    cellnum[$1]=$2;
                } else {
                    binid = ($1-xmin)/binsize;
                    if(binid>=0||binid<binnum){
                        t[int(binid)]+=$4;
                    }
                }
            }
            END {
              for(x=0;x<binnum;x++)
                  if(int(cellnum[x])==0) 
                      print 0;
                  else
                      print t[x]/cellnum[x]; 
            }'  ${indv}.cells_in_bin.txt $smes_file >$smes.txt
    else
        echo "failed for $smes_file"
    fi
done < $target_gene

touch tmp
awk -v binnum=$binnum '
    {}END{
        print "AP";
        for(x=0;x<binnum;x++)
            printf("AP_%d\n",x);
    }' tmp >column.txt

paste column.txt SMES*.txt >${indv}"_${binnum}.density.txt"
