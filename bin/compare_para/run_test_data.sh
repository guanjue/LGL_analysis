### set parameters
para_file1=01a2_S3V2_IDEAS_hg38_r1_rmh_rerun2.para.modified.para
para_file2=06a_S3V2_IDEAS_hg38_r3_withHg38Mm10prior.para.modified.para
num_of_feature=8
output_file=test.pdf

### run compare_two_para.R
Rscript compare_two_para.R $para_file1 $para_file2 $num_of_feature $output_file createGenomeTracks_np.R
