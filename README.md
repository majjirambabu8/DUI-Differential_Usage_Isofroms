# command 
Rscript dui_drimseq_v2.R \
  -g gencode.v33.primary_assembly.annotation.gtf \
  -s salmon_output_MUT_CNT_WT \ #salmon qunatiifcation output file
  -m metadata_Samples.tsv \
  -o Output_samples \
  -i sample_id \
  -c group \
  -r CNT



