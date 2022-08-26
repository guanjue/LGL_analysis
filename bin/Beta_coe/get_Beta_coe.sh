


### get cell-types in RNA-seq and IDEAS state matrices
declare -a common_ct=("B" "CLP" "CMP" "EOS" "ERY" "GMP" "LSK" "MEP" "MK" "MONc" "MONp" "MPP" "NEU" "NK" "CD4" "CD8")
### RNA-TPM matrix $RNA_TPM_file
declare -a rna_list=('LSK' 'LSK' 'ERY' 'ERY' 'CD4' 'CD4' 'CD8' 'CD8' 'B' 'B' 'CMP' 'CMP' 'MONp' 'MONp' 'NEU' 'NEU' 'MONc' 'MONc' 'GMP' 'GMP' 'CFUE' 'NK' 'NK' 'MK' 'MK' 'CLP' 'MPP' 'MPP' 'EOS' 'EOS' 'MEP' 'MEP' 'MK' 'MK' 'CLP' 'ERY' 'ERY' 'ERY' 'HUDEP1' 'HUDEP1' 'HUDEP2' 'HUDEP2' 'CD34' 'CD34')
### IDEAS matrix $IDEAS_state_file
declare -a sp_list=('AVE' 'B' 'B' 'CD34' 'CD34' 'CLP' 'CLP' 'CMP' 'CMP' 'EOS' 'EOS' 'ERY' 'ERY' 'GMP' 'GMP' 'HSC' 'HSC' 'HUDEP1' 'HUDEP1' 'HUDEP2' 'HUDEP2' 'K562' 'K562' 'LMPP' 'LMPP' 'MEP' 'MEP' 'MK' 'MK' 'MONc' 'MONc' 'MONp' 'MONp' 'MPP' 'MPP' 'NEU' 'NEU' 'NK' 'NK' 'CD4' 'CD4' 'CD8' 'CD8')
declare -a sp_ctrep_list=('AVE' 'B_B15_50' 'B_NC14_42' 'CD34_E_rep1' 'CD34_E_rep2' 'CLP_100266' 'CLP_100267' 'CMP_100246' 'CMP_100247' 'EOS_S006XEH2' 'EOS_S00BKK' 'ERY_S002R5' 'ERY_S002S3' 'GMP_100256' 'GMP_100257' 'HSC_100258' 'HSC_100259' 'HUDEP1_rep1' 'HUDEP1_rep2' 'HUDEP2_rep1' 'HUDEP2_rep2' 'K562_rep1' 'K562_rep2' 'LMPP_100268' 'LMPP_100269' 'MEP_Donor2596' 'MEP_Donor7256' 'MK_S004BTH2' 'MK_S00VHKH1' 'MONc_C0011IH1' 'MONc_C001UYH2' 'MONp_Prim_mon_C' 'MONp_Prim_mon_F' 'MPP_100272' 'MPP_100273' 'NEU_C0011IH2' 'NEU_C001UYH1' 'NK_S005YG' 'NK_S01E4WH0' 'T_CD4_S008H1' 'T_CD4_S009W4' 'T_CD8_C0066PH1' 'T_CD8_S00C2FH1')
### celltype list not used for Beta coefficients calculation
declare -a no_used_ct=('HUDEP1' 'HUDEP2' 'CD34')

Output_name_Beta_coefficient_mat='statep_rna_coe_heatmap.human.all.ccre.withcorfilter'
RNA_TPM_file_start='HumanVISION_RNAseq_hg38_genes_tpm'
All_RNA_gene_file_start='HumanVISION_RNAseq_hg38_gene'
cCRE_bed_file_start='S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU'
Working_folder='/Users/guanjuexiang/Documents/projects/analysis/04_05_2022_coe/coe_analysis_test'
rna_list_sample_num="${#rna_list[@]}"
sp_list_sample_num="${#sp_list[@]}"
no_used_ct_num="${#no_used_ct[@]}"

### state_rank for plotting
declare -a state_rank=(2 1 4 3 6 10 11 5 9 13 8 19 12 25 14 23 22 20 21 17 18 7 16 15 24)


state_num=25

statep_rna_coe_heatmap.human.all.ccre.withcorfilter.txt
HumanVISION_RNAseq_hg38_genes_tpm.idsort.protein_coding.txt
HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S[0:24].bed
HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S[0:24].bed
S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.S[0:24].mat.txt
HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withccreid.bed


for ct in "${common_ct[@]}"
do
echo "$ct"
time Rscript get_state_Beta_coefficients.human.R \
$ct \
$Output_name_Beta_coefficient_mat'.txt' \
$RNA_TPM_file_start'.idsort.protein_coding.txt' \
$All_RNA_gene_file_start'.idsort.protein_coding.Nkbupdownexp.S' \
$All_RNA_gene_file_start'.idsort.protein_coding.NHkbupdownexp.S' \
$cCRE_bed_file_start'.withid.S' \
$All_RNA_gene_file_start'.idsort.protein_coding.NHkbupdownexp.withccreid.bed' \
$Working_folder \
$rna_list_sample_num \
$sp_list_sample_num \
$no_used_ct_num \
$state_num \
"${rna_list[@]}" \
"${sp_list[@]}" \
"${no_used_ct[@]}" \
"${state_rank[@]}"
done

### get Beta coefficients for Human by taking ave beta cross all leave-one-out runs
Rscript get_ave.R \
'coe_score_no' $Output_name_Beta_coefficient_mat'.txt' $Output_name_Beta_coefficient_mat'.AVE' "${common_ct[@]}"





