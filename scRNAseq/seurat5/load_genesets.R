library(openxlsx)

geneset_path = "/data/project/genesets/"

# Neftel et al. Cell 2019
AC = read.table(paste0(geneset_path, "Neftel_AC.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
OPC = read.table(paste0(geneset_path, "Neftel_OPC.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
NPC1 = read.table(paste0(geneset_path, "Neftel_NPC1.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
NPC2 = read.table(paste0(geneset_path, "Neftel_NPC2.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
MES1 = read.table(paste0(geneset_path, "Neftel_MES1.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
MES2 = read.table(paste0(geneset_path, "Neftel_MES2.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]

neftel_list = list(AC, OPC, NPC1, NPC2, MES1, MES2)
names(neftel_list) = c("AC", "OPC", "NPC1", "NPC2", "MES1", "MES2")

# Sojka et al. Nature Cell Biology 2025
fetal_astrocyte = read.table(paste0(geneset_path, "Sojka_fetal_astrocyte.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
adult_astrocyte = read.table(paste0(geneset_path, "Sojka_adult_astrocyte.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
oligodendrocyte = read.table(paste0(geneset_path, "Sojka_oligodendrocyte.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
neuron = read.table(paste0(geneset_path, "Sojka_neuron.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
microglia = read.table(paste0(geneset_path, "Sojka_microglia.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]
endothelial_cells = read.table(paste0(geneset_path, "Sojka_endothelial_cells.txt"), sep="\t", header=T, stringsAsFactors=F)[[1]]

sojka_celltype_list = list(fetal_astrocyte, adult_astrocyte, oligodendrocyte, neuron, microglia, endothelial_cells)
names(sojka_celltype_list) = c("fetal_astrocyte", "adult_astrocyte", "oligodendrocyte", "neuron", "microglia", "endothelial_cells")

early = read.table(paste0(geneset_path, "Sojka_early100.txt"), sep="\t", header=F, stringsAsFactors=F)[[1]]
early_middle = read.table(paste0(geneset_path, "Sojka_early_middle100.txt"), sep="\t", header=F, stringsAsFactors=F)[[1]]
middle = read.table(paste0(geneset_path, "Sojka_middle100.txt"), sep="\t", header=F, stringsAsFactors=F)[[1]]
middle_late = read.table(paste0(geneset_path, "Sojka_middle_late100.txt"), sep="\t", header=F, stringsAsFactors=F)[[1]]
late = read.table(paste0(geneset_path, "Sojka_late100.txt"), sep="\t", header=F, stringsAsFactors=F)[[1]]
middle_coding = read.table(paste0(geneset_path, "sojka_middle_foundinscs_proteincoding.txt"), sep="\t", header=F, stringsAsFactors=F)[[1]]
early_middle_nonprol = read.table(paste0(geneset_path, "sojka_early_middle_non_proliferating_non_dnarepair.txt"), sep="\t", header=F, stringsAsFactors=F)[[1]]
middle_noncoding = setdiff(middle, middle_coding)
early_middle_prolif = setdiff(early_middle, early_middle_nonprol)

sojka_list = list(early=early, early_middle_nonprol=early_middle_nonprol, early_middle_prolif=early_middle_prolif, 
                middle_coding=middle_coding, middle_noncoding=middle_noncoding, middle_late=middle_late, late=late)

# Nomura et al. Nature Genetics 2025 
mp_RP = read.table(paste0(geneset_path, "metaprogram_RP.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_OPC = read.table(paste0(geneset_path, "metaprogram_OPC.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_CC = read.table(paste0(geneset_path, "metaprogram_CC.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_AC = read.table(paste0(geneset_path, "metaprogram_AC.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_Hypoxia = read.table(paste0(geneset_path, "metaprogram_Hypoxia.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_MES = read.table(paste0(geneset_path, "metaprogram_MES.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_NPC = read.table(paste0(geneset_path, "metaprogram_NPC.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_GPC = read.table(paste0(geneset_path, "metaprogram_GPC.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_ExN = read.table(paste0(geneset_path, "metaprogram_ExN.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_Stress1 = read.table(paste0(geneset_path, "metaprogram_Stress1.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
#mp_MIC = read.table(paste0(geneset_path, "metaprogram_MIC.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
#mp_LQ = read.table(paste0(geneset_path, "metaprogram_LQ.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_Cilia = read.table(paste0(geneset_path, "metaprogram_Cilia.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_NRGN = read.table(paste0(geneset_path, "metaprogram_NRGN.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
mp_Stress2 = read.table(paste0(geneset_path, "metaprogram_Stress2.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]

mp_list = list(mp_RP=mp_RP, mp_OPC=mp_OPC, mp_CC=mp_CC, mp_AC=mp_AC, mp_Hypoxia=mp_Hypoxia, mp_MES=mp_MES, mp_NPC=mp_NPC, 
            mp_GPC=mp_GPC, mp_ExN=mp_ExN, mp_Stress1=mp_Stress1, mp_Cilia=mp_Cilia, mp_NRGN=mp_NRGN, mp_Stress2=mp_Stress2)

# Hamed Nature 2025 Neural crest cell markers:
ncc = c("FOXD3", "ERBB3", "PLP1", "NGFR", "SOX10", "ETS1", "SPARC", "HES1")

# Read all Richards et al. Nature Cancer 2021 gene sets
richards_sets = read.xlsx(paste0(geneset_path, "Richards_2021_genesets.xlsx"), sheet=1)
richards_list = as.list(richards_sets)
richards_list = lapply(richards_list, function(x) x[!is.na(x)])

developmental_gsc = read.table(paste0(geneset_path, "Richards_2021_developmental_GSC.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
injury_gsc = read.table(paste0(geneset_path, "Richards_2021_injuryresponse_GSC.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
developmental_gsc_top = read.table(paste0(geneset_path, "Richards_2021_developmental_GSC_top.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]
injury_gsc_top = read.table(paste0(geneset_path, "Richards_2021_injuryresponse_GSC_top.txt"), header=F, sep="\t", stringsAsFactors=F)[[1]]

# Puchalski et al. Science 2018 
ivy = read.xlsx(paste0(geneset_path, "ivy_gap.xlsx"), sheet=1)

# Rautajoki et al. Acta Neuropathologica Communications 2023
hypoxia=c("CA9", "VEGFA", "ADM", "PDK1")
