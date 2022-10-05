library(chromoMap)
library(phytools)
library(dplyr)
library(ggplot2)
library(circlize)

#First step is to create the homology guide
#Read in the the orthogroup IDs, the Salar Gene to Orthogroup table and the Salar annnotations and merge
Ortho_names <- read.csv("Orthogroups.csv")
Ortho_data1 <- read.csv("Orthogroups.GeneCount.csv",header = TRUE)                    
Ortho_metadata <- inner_join(Ortho_data1,Ortho_names,by="Orthogroup")
Salar_gene_to_orthogroup <- read.csv("Salar_gene_to_orthogroup.csv")
Salar_annotations <- read.csv("Salar_annotation_EggNog.csv")
Salar_OG_annotations <- inner_join(Salar_gene_to_orthogroup,Salar_annotations, by="Salar")
Homology_guide <- inner_join(Ortho_metadata,Salar_OG_annotations,by="Orthogroup")
write.csv(Homology_guide, file = "Homology_guide.csv")

##############Ohnolog and Rediploidization IDs #####################

#Chinook
sum(Ortho_data1$Otsch)
#38425
Tsch_ohno <- filter(Ortho_data1, Otsch == 2* Esox & Otsch >= 1)
#9335
sum(Tsch_ohno$Otsch)
#19148
Tsch_redipping <- filter(Ortho_data1, Otsch < 2* Esox & Otsch > Esox & Esox >=1)
sum(Tsch_redipping$Otsch)
#1215
Tsch_expanding <- filter(Ortho_data1, Otsch > 2* Esox & Esox >= 1)
sum(Tsch_expanding$Otsch)
#1927
Tsch_contracting <- filter(Ortho_data1, Otsch > 0 & Otsch < Esox)
sum(Tsch_contracting$Otsch)
#661
Tsch_redip <- filter(Ortho_data1, Otsch == Esox & Otsch >= 1)
sum(Tsch_redip$Otsch)
#10157

##Coho 
sum(Ortho_data1$Okisu)
#40014
##Kisutch
Okisu_ohno <- filter(Ortho_data1, Okisu == 2* Esox & Okisu >= 1)
#10035
sum(Okisu_ohno$Okisu)
#20640
Okisu_rediping <- filter(Ortho_data1, Okisu < 2* Esox & Okisu >= Esox & Esox >=1)
sum(Okisu_rediping$Okisu)
#10542
Okis_expanding <- filter(Ortho_data1, Okisu > 2* Esox & Esox >=1)
sum(Okis_expanding$Okisu)
#2421
Okisu_redip <- filter(Ortho_data1, Okisu == Esox & Okisu >= 1)
#8933
sum(Okisu_redip$Okisu)
#9339
Okisu_contracting <- filter(Ortho_data1, Okisu > 0 & Okisu < Esox)
sum(Okisu_contracting$Okisu)
#638

##Sockeye
sum(Ortho_data1$Onerka)
#36503
Oner_ohno <- filter(Ortho_data1, Onerka == 2* Esox & Onerka >= 1)
#8466
#write.csv(Oner_ohno, file = "Onerka_ohnos.csv")
sum(Oner_ohno$Onerka)
#17266
Nerka_redipping <- filter(Ortho_data1, Onerka < 2* Esox & Onerka >= Esox & Esox >=1)
#write.csv(Nerka_other,file = "Nerka_rediploidizing.csv")
sum(Nerka_redipping$Onerka)
#10822
Nerka_expanding <- filter(Ortho_data1, Onerka > 2* Esox & Esox >= 1)
sum(Nerka_expanding$Onerka)
#1679
Nerka_contracting <- filter(Ortho_data1, Onerka > 0 & Onerka < Esox)
sum(Nerka_contracting$Onerka)
#764
Nerka_redip <- filter(Ortho_data1, Onerka == Esox & Onerka >=1 )
#write.csv(Nerka_redip, file = 'Nerka_redip.csv')
sum(Nerka_redip$Onerka)
#9969

##Chum
sum(Ortho_data1$Oketa)
#26153
Oket_ohno <- filter(Ortho_data1, Oketa == 2* Esox & Oketa >= 1)
#5773
sum(Oket_ohno$Oketa)
#11698
Keta_redipping <- filter(Ortho_data1, Oketa < 2* Esox & Oketa >= Esox & Esox >=1)
sum(Keta_redipping$Oketa)
#9690
Keta_expanding <- filter(Ortho_data1, Oketa > 2* Esox, Esox >= 1)
sum(Keta_expanding$Oketa)
#566
Keta_contracting <- filter(Ortho_data1, Oketa > 0 & Oketa < Esox)
sum(Keta_contracting$Oketa)
#658
Keta_redip <- filter(Ortho_data1, Oketa == Esox & Oketa >= 1)
sum(Keta_redip$Oketa)
#9208

##Rainbow Trout
sum(Ortho_data1$Omykiss)
#40620
Omyk_ohno <- filter(Ortho_data1, Omykiss == 2* Esox & Omykiss >= 1)
#10594
#write.csv(Omyk_ohno, file = "Mykiss_Ohns.csv")
sum(Omyk_ohno$Omykiss)
#21864
Mykiss_redipping <- filter(Ortho_data1, Omykiss < 2* Esox & Omykiss >= Esox & Esox >=1)
#write.csv(Mykiss_other, file = "Mykiss_rediploidizing.csv")
sum(Mykiss_redipping$Omykiss)
#10212
Mykiss_expanding <- filter(Ortho_data1, Omykiss > 2* Esox, Esox >=1)
sum(Mykiss_expanding$Omykiss)
#2238
Mykiss_contracting <- filter(Ortho_data1, Omykiss > 0 & Omykiss < Esox)
sum(Mykiss_contracting$Omykiss)
#633
Mykiss_redip <- filter(Ortho_data1, Omykiss == Esox & Omykiss >=1)
#write.csv(Mykiss_redip, file = "Mykiss_redip.csv")
sum(Mykiss_redip$Omykiss)
#8990

#Pink
sum(Ortho_data1$Ogorb)
#39324
Gorb_ohno <- filter(Ortho_data1, Ogorb == 2* Esox & Ogorb >= 1)
#9028
sum(Gorb_ohno$Ogorb)
#18540
Gorb_redipping <- filter(Ortho_data1, Ogorb < 2* Esox & Ogorb >= Esox & Esox >=1)
sum(Gorb_redipping$Ogorb)
#11104
Gorb_expanding <- filter(Ortho_data1, Ogorb > 2* Esox, Esox >=1)
sum(Gorb_expanding$Ogorb)
#3220
Gorb_contracting <- filter(Ortho_data1, Ogorb > 0 & Ogorb < Esox)
sum(Gorb_contracting$Ogorb)
#652
Gorb_redip <- filter(Ortho_data1, Ogorb == Esox & Ogorb >=1)
sum(Gorb_redip$Ogorb)
#10012

##Stats across Oncorhynchus
Ohnologs <- rbind(Okisu_ohno, Tsch_ohno, Oner_ohno, Oket_ohno,Omyk_ohno,Gorb_ohno)
Ohnologs_1 <- distinct(Ohnologs)
#12579
total_redip <- rbind(Okisu_redip,Tsch_redip,Nerka_redip,Keta_redip,Mykiss_redip,Gorb_redip)
total_redip1 <- distinct(total_redip)
#14776


##Create Figure 1
Ohnolog_barchart <- read.csv("Salmon_gene_classs_barplot.csv")
ggplot(Ohnolog_barchart, aes(fill=Class, y=Number, x=Species)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_viridis(discrete = T,direction =-1) +
  theme_classic()

#Find Overlap in retained duplicates for species tree in figure 2
Okisu_only_ohno <- filter(Ohnologs_1, Okisu == 2*Esox & Okisu >= 1 & Okisu != Oketa & Okisu != Omykiss & Okisu !=Onerka & Okisu != Otsch & Okisu!= Ogorb)
#214
Tsch_only_ohno <- filter(Ohnologs_1, Otsch == 2*Esox & Otsch>=1 & Otsch != Oketa & Otsch != Omykiss & Otsch !=Onerka & Otsch != Okisu & Otsch != Ogorb)
#233
Nerka_only_ohno <- filter(Ohnologs_1, Onerka == 2*Esox & Onerka >=1 & Onerka != Okisu & Onerka != Otsch & Onerka != Oketa & Onerka != Omykiss & Onerka != Ogorb)
#248
Keta_only_ohno <- filter(Ohnologs_1, Oketa == 2*Esox & Oketa >=1 & Oketa != Okisu & Oketa != Otsch & Oketa != Onerka & Oketa != Omykiss & Oketa != Ogorb)
#150
Mykisss_only_ohno <- filter(Ohnologs_1, Omykiss == 2*Esox & Omykiss >= 1 & Omykiss != Okisu & Omykiss != Otsch & Omykiss != Oketa & Omykiss != Onerka)
#582
Ogorb_only_ohno <- filter(Ohnologs_1, Ogorb == 2*Esox & Ogorb >= 1 & Ogorb != Okisu & Ogorb != Otsch & Ogorb != Oketa & Ogorb != Onerka & Ogorb != Omykiss)
#265

Keta_gorb_ohno <- filter(Ohnologs_1, Oketa == Ogorb & Oketa != Omykiss & Oketa != Esox & Oketa != Otsch & Oketa != Okisu & Oketa != Salar)
Keta_gorb_nerka <- filter(Keta_gorb_ohno, Oketa == Onerka)
Tsch_Okis <- filter(Ohnologs_1, Otsch == Okisu & Otsch != Onerka & Otsch != Oketa & Otsch !=Omykiss & Otsch != Salar & Otsch != Ogorb)
Tsch_Okis_Keta_gorb_nerka <- filter(Ohnologs_1, Otsch == Okisu & Otsch == Onerka & Otsch == Oketa & Otsch == Ogorb & Otsch !=Omykiss & Otsch != Salar)
Oncorhynchus_only_ohnos <- filter(Ohnologs_1, Otsch == Okisu & Otsch == Onerka & Otsch == Oketa & Otsch != Ogorb & Otsch ==Omykiss & Otsch != Salar)
All_ohno <- filter(Ohnologs_1,Otsch == 2*Esox & Otsch >= 1 & Otsch == Onerka & Otsch == Oketa & Otsch == Okisu & Otsch == Omykiss & Otsch == Ogorb)
#4164
All_Salar_ohno <- filter(All_ohno, Otsch == Salar)
#3741

#Find Overlap in rediploidized genes for species tree in figure 2
Okisu_only_redip <- filter(total_redip1, Okisu == Esox & Okisu != Oketa & Okisu != Omykiss & Okisu !=Onerka & Okisu != Otsch & Okisu != Ogorb)
#240
Tsch_only_redip <- filter(total_redip1, Otsch == Esox & Otsch != Oketa & Otsch != Omykiss & Otsch !=Onerka & Otsch != Okisu & Otsch != Ogorb)
#412
Nerka_only_redip <- filter(total_redip1, Onerka == Esox & Onerka != Okisu & Onerka != Otsch & Onerka != Oketa & Onerka != Omykiss & Onerka != Ogorb)
#565
Keta_only_redip <- filter(total_redip1, Oketa == Esox & Oketa != Okisu & Oketa != Otsch & Oketa != Onerka & Oketa != Omykiss & Oketa != Ogorb)
#1789
Mykisss_only_redip <- filter(total_redip1, Omykiss == Esox & Omykiss != Okisu & Omykiss != Otsch & Omykiss != Oketa & Omykiss != Onerka & Omykiss != Ogorb)
#222
Ogorb_only_redip <- filter(total_redip1, Ogorb == Esox & Ogorb != Okisu & Ogorb != Otsch & Ogorb !=Onerka & Ogorb != Oketa & Ogorb != Omykiss)
#487

Keta_gorb_redip <- filter(total_redip1, Oketa == Ogorb & Oketa != Omykiss & Oketa != Esox & Oketa != Otsch & Oketa != Okisu & Oketa != Salar)
Keta_gorb_nerka <- filter(Keta_gorb_redip, Oketa == Onerka)
Tsch_Okis_redip <- filter(total_redip1, Otsch == Okisu & Otsch != Onerka & Otsch != Oketa & Otsch !=Omykiss & Otsch != Salar & Otsch != Ogorb)
Redip_Tsch_Okis_Keta_gorb_nerka <- filter(total_redip1, Otsch == Okisu & Otsch == Onerka & Otsch == Oketa & Otsch == Ogorb & Otsch !=Omykiss & Otsch != Salar)
Oncorhynchus_only_redips <- filter(total_redip1, Otsch == Okisu & Otsch == Onerka & Otsch == Oketa & Otsch != Ogorb & Otsch ==Omykiss & Otsch != Salar)
All_redip <- filter(total_redip1,Otsch == Esox & Otsch == Onerka & Otsch == Oketa & Otsch == Okisu & Otsch == Omykiss & Otsch == Ogorb)
#5010
All_Salar <- filter(All_redip,Otsch == Esox & Otsch == Salar)
#4524

#Plot Tree for Figure 2
Salmonid_tree <- read.tree("SpeciesTree_rooted_node_labels.txt")
plot(Salmonid_tree)
#To generate the GO enrichments use the provided GO_MWU script
#Plot GO terms for figure2
mat <- as.data.frame(read.csv("Combined_chord_plot.csv"))
chordDiagram(mat)
#Fix up colors and labels for font in Illustrator





                             
