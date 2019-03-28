ggplot2::
  ggthemes::
  scales::
  ggrepel::
  RColorBrewer::
  grid::
  gridExtra:: 
  lattice::
  maftools::


## Read the files as MAF files 

Tumor_Germline <- read.delim("/Volumes/G-DRIVE mobile/6_samples_tumorgermline.maf.txt", header=TRUE, comment.char="#")
View(Tumor_Germline)
Tumor_Pipeline <- read.delim("/Volumes/G-DRIVE mobile/6_samples_tumoronly_pipeline_filtered.maf.txt", comment.char="#")
View(Tumor_Pipeline)
##For samples using germline filter
## Subsetting for the FFPE samples 
Tumor_Germline_subset <- Tumor_Germline[ which(Tumor_Germline$Tumor_Sample_Barcode == "CRC-2-FF"),]
View(Tumor_Germline_subset)
Tumor_Germline_subset1 <- Tumor_Germline[ which(Tumor_Germline$Tumor_Sample_Barcode == "CRC-1-FF"),]
View(Tumor_Germline_subset1)
Tumor_Germline_subset2 <- Tumor_Germline[ which(Tumor_Germline$Tumor_Sample_Barcode == "BreastCA-FF"),]
View(Tumor_Germline_subset2)
#combine all the FFPE samples 
Tumor_Germline_FF <- rbind(Tumor_Germline_subset, Tumor_Germline_subset1, Tumor_Germline_subset2)
View(Tumor_Germline_FF)
#Convert the data frame to MAF format
Tumor_Germline_FF_maf = maftools::read.maf(maf=Tumor_Germline_FF)

##For samples using pipeline as filter
##Subsetting for FFPE samples

Tumor_Pipeline_subset <- Tumor_Pipeline[ which(Tumor_Pipeline$Tumor_Sample_Barcode == "CRC-2-FF"),]
View(Tumor_Pipeline_subset)
Tumor_Pipeline_subset1 <- Tumor_Pipeline[ which(Tumor_Pipeline$Tumor_Sample_Barcode == "CRC-1-FF"),]
View(Tumor_Pipeline_subset1)
Tumor_Pipeline_subset2 <- Tumor_Pipeline[ which(Tumor_Pipeline$Tumor_Sample_Barcode == "BreastCA-FF"),]
View(Tumor_Pipeline_subset2)
##Combine all the FFPE samples
Tumor_Pipeline_FF <- rbind(Tumor_Pipeline_subset, Tumor_Pipeline_subset1, Tumor_Pipeline_subset2)
View(Tumor_Pipeline_FF)
##convert the dataframe to MAF format
Tumor_Pipeline_FF_maf = maftools::read.maf(maf=Tumor_Pipeline_FF)
#Comparison on Tumor_germline_maf and Tumor_Pipeline_maf
mafSummary(Tumor_Germline_FF_maf)
plotmafSummary(Tumor_Germline_FF_maf)
plotmafSummary(Tumor_Pipeline_FF_maf)
Tumor_Pipeline_FF_maf.titv = titv(maf = Tumor_Pipeline_FF_maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = Tumor_Pipeline_FF_maf.titv)




Tumor_Germline_maf <- read.maf(maf = "/Volumes/G-DRIVE mobile/6_samples_tumorgermline.maf.txt")
mafSummary(Tumor_Germline_maf)
Tumor_Pipeline_maf <- read.maf(maf = "/Volumes/G-DRIVE mobile/6_samples_tumoronly_pipeline_filtered.maf.txt")
mafSummary(Tumor_Pipeline_maf)
TMD_FFPE <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Chapter 4/https:/github.com/aditibagchi/Chapter_4.git/TMD_FFPE.csv")
TMD_FFPE <- as.data.frame(TMD_FFPE)
View(TMD_FFPE)
ggplot(data=TMD_FFPE, aes(x=TMD_FFPE$Case, y=TMD_FFPE$TMD, fill=Filtering)) +
  geom_bar(stat="identity", color = "black", position=position_dodge())+ 
  scale_y_continuous(name="Mutations per Mb", limits=c(0, 25)) +xlab("Cases")+
  scale_fill_brewer(palette = "Set1") + theme(axis.text.x = element_text(color = "black", size = 10, angle = 90)) + 
  geom_text(aes(label=TMD_FFPE$TMD), vjust=1.2, color="white",
            position = position_dodge(0.9), size=3.5)+
  theme_minimal() + ggtitle( "Source: FFPE" )
TMD_FROZEN <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/https:/github.com/aditibagchi/Chapter_4.git/TMD_FROZEN.csv")
ggplot(data=TMD_FROZEN, aes(x=TMD_FROZEN$Case, y=TMD_FROZEN$TMD, fill=Filtering_Strategy)) +
  geom_bar(stat="identity", color = "black", position=position_dodge())+ 
  scale_y_continuous(name="Mutations per Mb", limits=c(0, 25)) +xlab("Cases")+
  scale_fill_brewer(palette = "Set1") + theme(axis.text.x = element_text(color = "black", size = 10, angle = 90)) + 
  geom_text(aes(label=TMD_FROZEN$TMD), vjust=1.2, color="white",
            position = position_dodge(0.9), size=3.5)+
  theme_minimal() + ggtitle( "Source: FROZEN" )


##Mutation signatures
##Germline
#subset the data into "chr", "pos", "ref", "alt", "sample"
Tumor_Germline_FF_Sigs <- Tumor_Germline_FF[, c("Start_Position", "Chromosome", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")]
View(Tumor_Germline_FF_Sigs)
setnames(Tumor_Germline_FF_Sigs, old=c("Start_Position", "Chromosome", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"), new=c("pos","chr", "ref", "alt", "sample"))
str(Tumor_Germline_FF_Sigs)
Tumor_Germline_FF_Sigs = as.data.frame(Tumor_Germline_FF_Sigs)
sigs.input.Germline.FF = mut.to.sigs.input(mut.ref = Tumor_Germline_FF_Sigs, sample.id = "sample",chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = BSgenome.Hsapiens.UCSC.hg38)
View(sigs.input.Germline.FF)
BreastCA_FF  <- whichSignatures(tumor.ref = sigs.input.Germline.FF, sample.id = "BreastCA-FF",
                         signatures.ref = signatures.nature2013, associated = c(),
                         signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                         tri.counts.method = "default")
plotSignatures(BreastCA_FF)

##Mutationsignatures
##Filteringpipeline
##subset the data into "chr", "pos", "ref", "alt", "sample"

Tumor_Pipeline_FF_Sigs <- Tumor_Pipeline_FF[, c("Start_Position", "Chromosome", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")]
View(Tumor_Pipeline_FF_Sigs)
setnames(Tumor_Pipeline_FF_Sigs, old=c("Start_Position", "Chromosome", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"), new=c("pos","chr", "ref", "alt", "sample"))
sigs.input.Pipeline.FF = mut.to.sigs.input(mut.ref = Tumor_Pipeline_FF_Sigs, sample.id = "sample",chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = BSgenome.Hsapiens.UCSC.hg38)
View(sigs.input.Pipeline.FF)
CRC_2_FF  <- whichSignatures(tumor.ref = sigs.input.Pipeline.FF, sample.id = "CRC-2-FF",
                                signatures.ref = signatures.nature2013, associated = c(),
                                signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
                                tri.counts.method = "default")
plotSignatures(CRC_1_FF)
makePie(CRC_2_FF)

write.csv(sigs.input.Germline.FF, "sigs_input_germline_FFPE.csv")
write.csv(sigs.input.Pipeline.FF, "sigs_input_pipeline_FFPE.csv")  

## comparing the subsitution fractions 
subsitution_fractions <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Chapter 4/Chapter_4/subsitution_fractions.csv")
View(subsitution_fractions)
cols <- c("C>A" = "blue4", "C>G" = "grey66", "C>T" = "green4", "T>A" = "goldenrod2", "T>C" = "gray33", "T>G" = "red3")

X = ggplot(data= subsitution_fractions, aes(x=subsitution_fractions$Filtering, y=subsitution_fractions$Fraction, fill=Subsitution)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Subsitutions") + theme_minimal()+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90)) 
S= X + facet_wrap(~Tumor)


