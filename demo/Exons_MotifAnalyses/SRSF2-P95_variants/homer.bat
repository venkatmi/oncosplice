#BSUB -L /bin/bash
#BSUB -W 50:00
#BSUB -M 10000
#BSUB -n 1
#BSUB -e /users/ches2d/lsf_logs/%J.err
#BSUB -u Xiaoting.Chen@cchmc.org
# execute program
source ~/.bashrc
module load gcc/4.8.1

cd /data/weirauchlab/team/ches2d/OngoingProjects/Meenakshi_enrichment_analysis_Jan2018/Round1/Exons/SRSF2-P95

findMotifsGenome.pl input.bed hg19 output1/ -size given -rna -mknown /data/weirauchlab/web/basedata/HomerData/customizedRNAs/v0.60/human_cisbp060.motif
findMotifsGenome.pl input.bed hg19 output2/ -size given -rna -mknown /data/weirauchlab/web/basedata/HomerData/customizedRNAs/v0.60/human_cisbp060.motif
findMotifsGenome.pl input.bed hg19 output3/ -size given -rna -mknown /data/weirauchlab/web/basedata/HomerData/customizedRNAs/v0.60/human_cisbp060.motif
findMotifsGenome.pl input.bed hg19 output4/ -size given -rna -mknown /data/weirauchlab/web/basedata/HomerData/customizedRNAs/v0.60/human_cisbp060.motif
findMotifsGenome.pl input.bed hg19 output5/ -size given -rna -mknown /data/weirauchlab/web/basedata/HomerData/customizedRNAs/v0.60/human_cisbp060.motif
#findMotifsGenome.pl input.bed hg19 output_TF/ -size given -nomotif -mknown /data/weirauchlab/web/basedata/HomerData/customizedRNAs/v0.60/human_cisbp060.motif
findMotifsGenome.pl input.bed hg19 output_TF_strand/ -size given -nomotif -norevopp -mknown    /data/weirauchlab/web/basedata/HomerData/customizedRNAs/v0.60/human_cisbp060.motif
