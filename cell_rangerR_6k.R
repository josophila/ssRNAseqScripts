#plotting genes onto the gransden 6k data


library(cellrangerRkit)

cellranger_pipestance_path = "/Users/Josophila/Desktop/gransden_6k/"
gbm = load_cellranger_matrix(cellranger_pipestance_path)
analysis_results = load_cellranger_analysis_results(cellranger_pipestance_path)


#look at tSNE
tsne_proj = analysis_results$tsne
visualize_umi_counts(gbm,tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(3,4), marker_size = 0.5)

#Filter unexpressed genes
use_genes = get_nonzero_genes(gbm)

#normalize the UMI counts for each barcode
gbm_bcnorm = normalize_barcode_sums_to_median(gbm[use_genes,])

#This is the log-transformed, normalized gene barcode matrix
gbm_log = log_gene_bc_matrix(gbm_bcnorm, base = 10)

#After - number of genes detected and number of cells detected
print(dim(gbm_log))

#Visualizing signatures of known gene markers
#Import synonym dataframe
rawsyn = read.csv('C:/Users/Josophila/Desktop/00_PhD/11_bioinformatics_tools/000_moss_genome/Ppatens/annotation/PhytozomeV11_download/annotation/syn_simple_prog.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
colnames(rawsyn) = c("version3_ID", "version1.6_ID")

#Mike asked for a few gene names.  Pp1s171_161V6 is a mis-annotation.
mjscaul = 'Pp1s171_161V6|Pp1s338_11V6'
mjs_tip = 'Pp1s346_19V6|Pp1s28_242V6|Pp1s103_27V6'

caulconv = grep(rawsyn$version1.6_ID, pattern = mjscaul)
tipconv = grep(rawsyn$version1.6_ID, pattern = mjs_tip)

mjscaulv3 = rawsyn[caulconv,]$version3_ID
mjs_tipv3 = rawsyn[tipconv,]$version3_ID
mjsgenes = c(mjscaulv3, mjs_tipv3)

#visualize expression of mike's chosen genes
tsne_proj = analysis_results$tsne


mhf_top_tip = "Pp1s247_15V6|Pp1s385_3V6|Pp1s218_28V6|Pp1s215_82V6|Pp1s9_179V6|Pp1s169_74V6|Pp1s242_95V6|Pp1s225_73V6|Pp1s31_253V6"
mhfconv = grep(rawsyn$version1.6_ID, pattern = mhf_top_tip)
mhftips = rawsyn[mhfconv,]$version3_ID



mhf_high_tip = "Pp1s249_69V6|Pp1s158_100V6|Pp1s249_16V6|Pp1s82_6V6|Pp1s472_6V6|Pp1s184_89V6|Pp1s237_75V6|Pp1s54_29V6|Pp1s191_56V6"
mhfentipconv = grep(rawsyn$version1.6_ID, pattern = mhf_high_tip)
mhf_en_tips = rawsyn[mhfentipconv,]$version3_ID


#from_GENE_Atlas_data - preferentially expressed genes from each
genat_chlor = 'Pp1s16_216V6|Pp1s108_49V6|Pp1s199_50V6|Pp1s70_55V6|Pp1s124_54V6|Pp1s187_117V6|Pp1s226_76V6|Pp1s124_56V6|Pp1s148_17V6|Pp1s181_103V6|Pp1s50_57V2|Pp1s1_34V6|Pp1s69_132V6|Pp1s15_193V6|Pp1s280_53V6|Pp1s139_78V6'


genat_rhiz= 'Pp1s115_69V6|Pp1s122_19V6|Pp1s147_139V6|Pp1s149_101V6|Pp1s15_359V6|Pp1s162_149V6|Pp1s167_11V6|Pp1s174_54V6|Pp1s19_11V6|Pp1s199_79V6|Pp1s2_334V6|Pp1s20_266V6|Pp1s202_94V6|Pp1s217_21V6|Pp1s24_14V6|Pp1s249_69V6'

genat_caul = 'Pp1s46_140V6|Pp1s3_547V6|Pp1s46_32V6|Pp1s96_59V6|Pp1s86_93V6|Pp1s5_311V6|Pp1s233_46V6|Pp1s926_1V6|Pp1s143_38V6|Pp1s159_139V6|Pp1s59_246V6|Pp1s121_60V6|Pp1s61_299V6|Pp1s27_207V6|Pp1s57_140V6|Pp1s100_89V6'

genat_chlorEn = 'Pp1s66_48V6|Pp1s142_72V6|Pp1s66_47V6|Pp1s243_73V6|Pp1s76_48V6|Pp1s51_303V6|Pp1s107_173V6|Pp1s55_253V6|Pp1s133_16V6|Pp1s375_54V6|Pp1s68_116V6|Pp1s82_3V6|Pp1s56_177V6|Pp1s211_112V6|Pp1s68_12V6|Pp1s2_770V6|Pp1s41_267V6|Pp1s82_6V6|Pp1s370_30V6|Pp1s56_13V6|Pp1s21_274V6|Pp1s217_66V6|Pp1s56_17V6|Pp1s465_25V6'

chlorENgrep = grep(rawsyn$version1.6_ID, pattern = genat_chlorEn)
chlorEnv3 = rawsyn[chlorENgrep,'version3_ID']
chlorEnv3
chlorEnv3 = chlorEnv3[-c(13,11,9,1)]

caulEn = 'Pp1s428_13V6|Pp1s47_88V6|Pp1s164_49V6|Pp1s378_15V6|Pp1s17_180V6|Pp1s14_2V6|Pp1s240_2V6|Pp1s183_29V6|Pp1s328_44V6|Pp1s23_51V6|Pp1s182_93V6|Pp1s27_232V6|Pp1s68_110V6|Pp1s309_69V6|Pp1s346_19V6|Pp1s120_113V2|Pp1s195_100V6|Pp1s47_244V6|Pp1s111_186V6|Pp1s82_83V6|Pp1s25_61V6|Pp1s143_12V6|Pp1s68_127V6|Pp1s122_57V6|Pp1s295_3V6|Pp1s333_41V6|Pp1s334_84V6|Pp1s309_84V6'
caulEngrep = grep(rawsyn$version1.6_ID, pattern = caulEn)
caulEnV3 = rawsyn[caulEngrep, 'version3_ID']

caulEnV3
caulEnV3 = caulEnV3[-c(23, 11, 6)]

visualize_gene_markers(gbm_log, caulEnV3[1:16], tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5), title = "Top Enriched in Caulonema")

genat_chlor_match = grep(rawsyn$version1.6_ID, pattern = genat_chlor)
genat_caul_match = grep(rawsyn$version1.6_ID, pattern = genat_caul)
genat_rhiz_match = grep(rawsyn$version1.6_ID, pattern = genat_rhiz)

genat_chlor_pref = rawsyn[genat_chlor_match,'version3_ID']
genat_caul_pref = rawsyn[genat_caul_match, 'version3_ID']
genat_rhiz_pref = rawsyn[genat_rhiz_match,'version3_ID']

genat_rhiz_pref
genat_rhiz_pref = genat_rhiz_pref[-c(9, 7)]

genat_caul_pref
genat_caul_pref = genat_caul_pref[-c(13,11,4)]
 


###Upregulated in cold 4hr response
up_cold = 'Pp1s142_63V6|Pp1s142_72V6|Pp1s26_152V6|Pp1s12_415V6|Pp1s41_47V6|Pp1s92_84V6|Pp1s55_66V6|Pp1s142_76V6|Pp1s128_100V6|Pp1s194_80V6|Pp1s660_1V6|Pp1s78_27V6|Pp1s147_37V6|Pp1s102_107V6|Pp1s51_40V6|Pp1s156_153V6|Pp1s122_25V6|Pp1s100_144V6|Pp1s171_52V6'
up_coldgrep = grep(rawsyn$version1.6_ID, pattern = up_cold)
upcoldv3 = rawsyn[up_coldgrep,'version3_ID']

upcoldv3 = upcoldv3[-8]


up_aba4hr = 'Pp1s341_25V6|Pp1s17_269V6|Pp1s153_46V6|Pp1s46_232V6|Pp1s465_25V6|Pp1s232_45V6|Pp1s210_89V6|Pp1s246_98V6|Pp1s17_101V6|Pp1s419_7V6|Pp1s18_14V6|Pp1s34_430V6|Pp1s257_115V6|Pp1s35_60V6|Pp1s4_386V6|Pp1s141_116V6|Pp1s41_47V6|Pp1s93_48V6|Pp1s7_419V6|Pp1s54_108V6|Pp1s5_162V6'
upaba4_grep = grep(rawsyn$version1.6_ID, pattern = up_aba4hr)
upaba_4hr_v3 = rawsyn[upaba4_grep, 'version3_ID']

upaba_4hr_v3 = upaba_4hr_v3[-17]
pdf(file = './tsne_gene_plots_1-11-17.pdf')

visualize_gene_markers(gbm_log, mjscaulv3, tsne_proj[c("TSNE.1","TSNE.2")], limits =c(0,1.5), title = 'suggested caulonema gene')
visualize_gene_markers(gbm_log, mjs_tipv3, tsne_proj[c("TSNE.1","TSNE.2")], limits =c(0,1.5), title = 'suggested tip gene')

visualize_gene_markers(gbm_log, mhf_en_tips, tsne_proj[c("TSNE.1","TSNE.2")], limits = c(0,1.5), title = 'Highest enriched tip genes MHF')

visualize_gene_markers(gbm_log, mhftips, tsne_proj[c("TSNE.1","TSNE.2")], limits = c(0,1.5), title = 'Most significant Tip MHF')

visualize_gene_markers(gbm_log, genat_chlor_pref, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5), title = "Preferentially in Chloronema")

visualize_gene_markers(gbm_log, chlorEnv3, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5), title = "Top Enriched chloronema")

visualize_gene_markers(gbm_log, genat_caul_pref, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5), title = "Preferentially in Caulonema")

visualize_gene_markers(gbm_log, caulEnV3[1:16], tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5), title = "Top Enriched in Caulonema")

visualize_gene_markers(gbm_log, genat_rhiz_pref, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5), title = "Preferentially in Rhizoids")

visualize_gene_markers(gbm_log, upcoldv3, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5), title = "Upregulated in Cold after 4 hr")

visualize_gene_markers(gbm_log, upaba_4hr_v3, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5), title = "Upregulated in ABA after 4 hr")

dev.off()

