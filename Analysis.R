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
View(Tumor_Germline)
Tumor_Germline_maf <- read.maf(maf = "/Volumes/G-DRIVE mobile/6_samples_tumorgermline.maf.txt")
mafSummary(Tumor_Germline_maf)
Tumor_Pipeline_maf <- read.maf(maf = "/Volumes/G-DRIVE mobile/6_samples_tumoronly_pipeline_filtered.maf.txt")
mafSummary(Tumor_Pipeline_maf)
TMD_FFPE <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/https:/github.com/aditibagchi/Chapter_4.git/TMD_FFPE.csv")
TMD_FFPE <- as.data.frame(TMD_FFPE)
View(TMD_FFPE)
ggplot(data=TMD_FFPE, aes(x=TMD_FFPE$Case, y=TMD_FFPE$TMD, fill=Filtering_Strategy)) +
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

                                                                                                                                        