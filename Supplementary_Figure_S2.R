## Code to generate Supplementary Figure 2 from 
## "AACR Project GENIE: 100,000 cases and beyond", 2021, 
## AACR Project GENIE Consortium, Genomics and Analysis Working Group

library(reshape2)
library(plyr)
library(ggplot2)

data_clinical_sample <- read.table("data_clinical_sample_9.1-public.txt",
              sep = "\t", header = TRUE, comment.char = "",quote = "",
              skip = 4)

Tumour_type_key <- setNames(data_clinical_sample$CANCER_TYPE,
                            data_clinical_sample$SAMPLE_ID)

oncotree_key <- setNames(data_clinical_sample$ONCOTREE_CODE,
                            data_clinical_sample$SAMPLE_ID)

maf <- read.table("data_mutations_extended_9.1-public.txt",
                         sep = "\t", header = TRUE, 
                         comment.char = "", quote = "")

mutation_counts <- as.data.frame(table(maf$Tumor_Sample_Barcode))

colnames(mutation_counts) <- c("Sample","Mutation_Count")

## Add Cancer type and oncotree code to maf file for oncoKB annotator to 
## accurately annotate variants.
maf$CANCER_TYPE  <-  Tumour_type_key[maf$Tumor_Sample_Barcode]
maf$ONCOTREE_CODE <- oncotree_key[maf$Tumor_Sample_Barcode]

write.table(maf, "data_mutations_extended_9.1-public_with_tumour_type.txt",
           sep = "\t", row.names=F, quote = F)

## Run below at the command line to annotate maf with oncokb annotator:
## https://github.com/oncokb/oncokb-annotator

'
python3 MafAnnotator.py -i data_mutations_extended_9.1-public_with_tumour_type.txt \
-o data_mutations_extended_9.1-public_oncokb_annotated.maf \
-b $TOKEN
'

## Read in oncokb annotated maf
maf <- read.table("data_mutations_extended_9.1-public_oncokb_annotated.maf",
                  sep = "\t", header = TRUE, 
                  comment.char = "", quote = "")

## Summarize the number of putative "driver" vs. "passenger" mutations

maf$ONCOGENIC[maf$ONCOGENIC==""] <- "Unknown"

counts <- ddply(maf, .(maf$Tumor_Sample_Barcode, maf$ONCOGENIC), nrow)

colnames(counts) <- c("Tumor_Sample_Barcode", "Category", "Freq")

counts$Category <- gsub(" ","_",counts$Category)

samples <- unique(counts$Tumor_Sample_Barcode)

## Summarize the number of "Driver" mutations per sample.  For this purpose, 
## "Driver" is defined as having an "Oncogenic", "Likely Oncogenic", "Predicted Oncogenic" 
## or "Resistance" labels from oncoKB

counts$Category <- ifelse(counts$Category == "Oncogenic" | 
                            counts$Category == "Likely_Oncogenic" |
                            counts$Category == "Predicted_Oncogenic" |
                            counts$Category == "Resistance", "Driver", "Non_Driver")

cats <- unique(counts$Category)

mutation_counts <- matrix(data = 0,
                          nrow = length(samples),
                          ncol = length(cats),
                          dimnames = list(samples,cats))

mutation_counts[cbind(counts$Tumor_Sample_Barcode,counts$Category)] <- counts$Freq
mutation_counts <-as.data.frame(mutation_counts)
mutation_counts$Mutation_Count <- as.numeric(apply(mutation_counts,1,sum))

## Add non-mutated samples (these aren't included in the maf but are in the sample file)
No_mut_samples <- setdiff(data_clinical_sample$SAMPLE_ID, maf$Tumor_Sample_Barcode)
cats <- c(cats,"Mutation_Count")

No_mut_frame <- as.data.frame(matrix(data = 0,
                              nrow = length(No_mut_samples),
                              ncol = length(cats),
                              dimnames = list(No_mut_samples,cats)))

mutation_counts <- rbind(mutation_counts,
                         No_mut_frame)

mutation_counts$Cancer_Type <-  Tumour_type_key[row.names(mutation_counts)]

## Summarize frequency of samples with: 
## a) at least one driver, 
## b) only "non-drivers" or 
## c) no mutations in each cancer type

mutation_counts$Class <- ifelse(mutation_counts$Driver !=0, "With_Drivers",
                               ifelse(mutation_counts$Mutation_Count !=0,"Only_Non_Drivers",
                                      "No_Mutations"))
  
bar_data <- ddply(mutation_counts, .(mutation_counts$Class, mutation_counts$Cancer_Type), nrow)

colnames(bar_data) <- c("Class","Cancer_Type","Frequency")

## Plot the Top 30 cancer types with the highest sample count
cancer_type_summ <- as.data.frame(table(data_clinical_sample$CANCER_TYPE))
cancer_type_summ <- cancer_type_summ[order(cancer_type_summ$Freq,decreasing = T),]
colnames(cancer_type_summ) <- c("Cancer_Type", "Freq")
## remove "UNKNOWNS" Cancer Type
cancer_type_summ <- subset(cancer_type_summ, Cancer_Type !="UNKNOWN")
row.names(cancer_type_summ) <- cancer_type_summ$Cancer_Type

cancer_type_30 <- as.character(cancer_type_summ$Cancer_Type[1:30])

bar_data <- subset(bar_data,bar_data$Cancer_Type %in% cancer_type_30)

## Order based on proportion of "With Driver" samples 
cancer_order <- as.data.frame(table(mutation_counts$Cancer_Type[mutation_counts$Class=="With_Drivers"]))
colnames(cancer_order) <- c("Cancer_Type","With_Driver_Count")
cancer_order <- subset(cancer_order,cancer_order$Cancer_Type %in% cancer_type_30)

cancer_order <- merge(cancer_order,cancer_type_summ,by="Cancer_Type")
cancer_order <- cancer_order[order(cancer_order$With_Driver_Count/cancer_order$Freq),]

bar_data$Cancer_Type <- factor(bar_data$Cancer_Type, levels=cancer_order$Cancer_Type)

bar_data$Total <- cancer_type_summ[as.character(bar_data$Cancer_Type),"Freq"]

pdf(paste0("~/OneDrive - UHN/PughLab_Work/GENIE/Figures/",
           "Driver_mutation_by_Tumour_Type_v9.1.pdf", sep=""),
    height = 6, width=17)
ggplot(bar_data, 
       aes(fill=Class, y=Frequency, x= Cancer_Type )) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1.1,face="bold"),
        axis.text.y = element_text(face="bold"),
        axis.line = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank()) + 
  scale_y_continuous(labels = scales::percent) +
  annotate("text", x = 1:30, y = 1.03, label = cancer_order$Freq) +
  annotate("text", x = 31.5, y = 1.03, label = "(n samples)") +
  coord_cartesian(clip = "off") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(x="", y="Percentage")
dev.off()
