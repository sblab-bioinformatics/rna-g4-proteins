#!/bin/bash

path='.'

gene_sets_user_names=( 'ER' 'cellular_component' 'hallmarkGeneSets' 'c2_curatedGeneSets' 'oncogenicSignatures'  )
gene_sets_DB_names=( 'ER.gmt' 'c5.cc.v6.2.symbols.gmt' 'h.all.v6.2.symbols.gmt' 'c2.all.v6.2.symbols.gmt' 'c6.all.v6.2.symbols.gmt'  )


for i in `seq ${#gene_sets_user_names[@]}`
do
   j=$(( $i -1 ))
   gene_sets_user_name=${gene_sets_user_names[$j]}
   gene_sets_DB_name=${gene_sets_DB_names[$j]}
   echo "$gene_sets_user_name == $gene_sets_DB_name"
   gmx_file=$path/MSigDB/$gene_sets_DB_name



   for each_gene_set in $( ls ${path}/*.rnk)
   do
      echo $each_gene_set
      contrast=$(basename $each_gene_set | sed 's/_v2.0.rnk//'  )
      echo $contrast
      rankfile=$each_gene_set

      if [ ! -d "$path/output/$gene_sets_user_name" ]
      then
         mkdir -p $path/output/$gene_sets_user_name
      fi
      outputdir=$path/output/$gene_sets_user_name

      java -cp ~/gsea-3.0.jar -Xmx2G xtools.gsea.GseaPreranked \
           -gmx $gmx_file \
           -collapse false \
           -nperm 1000 \
           -rnk $rankfile \
           -scoring_scheme weighted \
           -rpt_label $contrast \
           -include_only_symbols true \
           -make_sets true -plot_top_x 40 \
           -rnd_seed 1 \
           -set_max 500 \
           -set_min 5 \
           -zip_report false \
           -out $outputdir \
           -gui false \
           -create_svgs true
   done
done
