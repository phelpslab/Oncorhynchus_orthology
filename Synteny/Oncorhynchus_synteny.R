##Chinook Synteny

Tsch_protins <- read.csv("Tsch_proteins.csv")
Orthogroup_IDs <- read.csv("Orthogroup_protein_IDs.csv")
Tsch_conserved_proteins <- inner_join(Orthogroup_IDs,Ortho_data1,by="Orthogroup")
write.csv(Tsch_conserved_proteins, file = "Tsch_conserved_protein_list.csv")


Tsch_annotations <- cbind(Tsch_test,Tsch_protein)
Tsch_proteins <- filter(Tsch_annotations, Tsch_test == TRUE)

Chrom1_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG01' & Start >= 8200030 & Start <= 37365193)
Chrom1_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG01' & Start >= 44689526 & Start <= 44793320)
Chrom1_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG01' & Start >= 51287382 & Start <= 74786573)
Chrom2_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG02' & Start >= 23446809 & Start <= 23591168)
Chrom2_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG02' & Start >= 4806009 & Start <= 45591258)
Chrom3_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG03' & Start >= 46067475 & Start <= 65244049)
Chrom3_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG03' & Start >= 6799445 & Start <= 22097843)
Chrom3_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG03' & Start >= 11756496 & Start <= 17181698)
Chrom3_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG03' & Start >= 44849857 & Start <= 44849857)
Chrom4_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG04' & Start >= 56728251 & Start <= 69310764)
Chrom4_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG04' & Start >= 2857796 & Start <= 69310764)
Chrom4_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG04' & Start >= 5194345 & Start <= 646289)
Chrom4_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG04' & Start >= 70595607 & Start <= 70753918)
Chrom4_block5 <- filter(Tsch_protins, Chrom == 'linkage group LG04' & Start >= 21571784 & Start <= 25059591)
Chrom5_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG05' & Start >= 22283161 & Start <= 25545682)
Chrom5_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG05' & Start >= 42394974 & Start <= 42848264)
Chrom5_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG05' & Start >= 48929649 & Start <= 72294367)
Chrom5_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG05' & Start >= 12495673 & Start <= 41371419)
Chrom6_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG06' & Start >= 6657351 & Start <= 6772118)
Chrom6_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG06' & Start >= 5250226 & Start <= 10654839)
Chrom6_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG06' & Start >= 53034602 & Start <= 71575415)
Chrom6_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG06' & Start >= 15134872 & Start <= 46237325)
Chrom7_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG07' & Start >= 27738763 & Start <= 27866041)
Chrom7_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG07' & Start >= 8606258 & Start <= 8736023)
Chrom7_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG07' & Start >= 33903579 & Start <= 79712201)
Chrom8_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG08' & Start >= 9149916 & Start <= 29514087)
Chrom8_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG08' & Start >= 46323767 & Start <= 60836874)
Chrom9_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG09' & Start >= 14806976 & Start <= 36968390)
Chrom9_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG09' & Start >= 65321468 & Start <= 67127944)
Chrom9_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG09' & Start >= 68698376 & Start <= 78021929)
Chrom9_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG09' & Start >= 82673855 & Start <= 82540255)
Chrom10_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG10' & Start >= 16754272 & Start <= 1874393)
Chrom10_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG10' & Start >= 2823595 & Start <= 28537602)
Chrom10_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG10' & Start >= 33756987 & Start <= 51372922)
Chrom12_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG12' & Start >= 20417484 & Start <= 44847904)
Chrom12_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG12' & Start >= 34145442 & Start <= 34257985)
Chrom12_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG12' & Start >= 39157814 & Start <= 39441080)
Chrom12_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG12' & Start >= 52324737 & Start <= 67483126)
Chrom12_block5 <- filter(Tsch_protins, Chrom == 'linkage group LG12' & Start >= 72749427 & Start <= 73089051)
Chrom13_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG13' & Start >= 14950143 & Start <= 15233409)
Chrom13_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG13' & Start >= 16296612 & Start <= 33797824)
Chrom13_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG13' & Start >= 37178227 & Start <= 41610013)
Chrom13_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG13' & Start >= 41785862 & Start <= 41900630)
Chrom13_block5 <- filter(Tsch_protins, Chrom == 'linkage group LG13' & Start >= 46770077 & Start <= 59638676)
Chrom14_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG14' & Start >= 21903405 & Start <= 39669713)
Chrom14_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG14' & Start >= 2534323 & Start <= 13280901)
Chrom15_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG15' & Start >= 6268372 & Start <= 11668493)
Chrom15_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG15' & Start >= 12450838 & Start <= 13720076)
Chrom15_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG15' & Start >= 15618658 & Start <= 15733705)
Chrom15_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG15' & Start >= 16828362 & Start <= 16971153)
Chrom15_block5 <- filter(Tsch_protins, Chrom == 'linkage group LG15' & Start >= 20589197 & Start <= 22063273)
Chrom15_block6 <- filter(Tsch_protins, Chrom == 'linkage group LG15' & Start >= 24916774 & Start <= 28771994)
Chrom16_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG16' & Start >= 34449793 & Start <= 50836574)
Chrom16_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG16' & Start >= 57271526 & Start <= 57693722)
Chrom16_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG16' & Start >= 13739096 & Start <= 19046089)
Chrom16_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG16' & Start >= 26339791 & Start <= 29347345)
Chrom16_block5 <- filter(Tsch_protins, Chrom == 'linkage group LG16' & Start >= 57271526 & Start <= 57693722)
Chrom17_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG17' & Start >= 8898903 & Start <= 19192947)
Chrom18_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG18' & Start >= 9485311 & Start <= 15704342)
Chrom18_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG18' & Start >= 22662135 & Start <= 22774151)
Chrom19_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG19' & Start >= 19373921 & Start <= 46245744)
Chrom20_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG20' & Start >= 22279729 & Start <= 22385373)
Chrom20_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG20' & Start >= 3960432 & Start <= 36910044)
Chrom21_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG21' & Start >= 3426122 & Start <= 4155862)
Chrom21_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG21' & Start >= 9960240 & Start <= 20172701)
Chrom21_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG21' & Start >= 27026740 & Start <= 27205202)
Chrom22_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG22' & Start >= 6871097 & Start <= 17896558)
Chrom22_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG22' & Start >= 28890786 & Start <= 28766866)
Chrom23_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG23' & Start >= 3387692 & Start <= 3529796)
Chrom23_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG23' & Start >= 8841094 & Start <= 9013673)
Chrom23_block3 <- filter(Tsch_protins, Chrom == 'linkage group LG23' & Start >= 10758765 & Start <= 10915862)
Chrom23_block4 <- filter(Tsch_protins, Chrom == 'linkage group LG23' & Start >= 11348562 & Start <= 19134363)
Chrom24_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG24' & Start >= 21896194 & Start <= 22006735)
Chrom24_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG24' & Start >= 9622013 & Start <= 19017057)
Chrom25_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG25' & Start >= 6046399 & Start <= 33015369)
Chrom26_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG26' & Start >= 549814 & Start <= 13515016)
Chrom26_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG26' & Start >= 26858855 & Start <= 35019694)
Chrom27_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG27' & Start >= 2322447 & Start <= 17539852)
Chrom28_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG28' & Start >= 6123273 & Start <= 19869005)
Chrom28_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG28' & Start >= 21108178 & Start <= 25056005)
Chrom29_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG29' & Start >= 15522487 & Start <= 25604190)
Chrom30_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG30' & Start >= 16920121 & Start <= 36132668)
Chrom31_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG31' & Start >= 52837511 & Start <= 78590339)
Chrom32_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG32' & Start >= 9345424 & Start <= 9473630)
Chrom33_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG33' & Start >= 13006417 & Start <= 22124850)
Chrom33_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG33' & Start >= 27142002 & Start <= 32818981)
Chrom34_block1 <- filter(Tsch_protins, Chrom == 'linkage group LG34' & Start >= 1719280 & Start <= 1893722)
Chrom34_block2 <- filter(Tsch_protins, Chrom == 'linkage group LG34' & Start >= 6596123 & Start <= 6703084)

Tsch_syntenic_genes <- rbind(Chrom1_block1 ,Chrom1_block2 ,Chrom1_block3 ,Chrom2_block1 ,Chrom2_block2 ,Chrom3_block1 ,Chrom3_block2 ,Chrom3_block3 ,Chrom3_block4 ,Chrom4_block1 ,Chrom4_block2 ,Chrom4_block3 ,Chrom4_block4 ,Chrom4_block5 ,Chrom5_block1 ,Chrom5_block2 ,Chrom5_block3 ,Chrom5_block4 ,Chrom6_block1 ,Chrom6_block2 ,Chrom6_block3 ,Chrom6_block3 ,Chrom7_block1 ,Chrom7_block2 ,Chrom7_block3 ,Chrom8_block1 ,Chrom8_block2 ,Chrom9_block1 ,Chrom9_block2 ,Chrom9_block3 ,Chrom9_block4 ,Chrom10_block1 ,Chrom10_block2 ,Chrom10_block3 ,Chrom12_block1 ,Chrom12_block2 ,Chrom12_block3 ,Chrom12_block4 ,Chrom12_block5 ,Chrom13_block1 ,Chrom13_block2 ,Chrom13_block3 ,Chrom13_block4 ,Chrom13_block5 ,Chrom14_block1 ,Chrom14_block2 ,Chrom15_block1 ,Chrom15_block2 ,Chrom15_block3 ,Chrom15_block4 ,Chrom15_block5 ,Chrom15_block6 ,Chrom16_block1 ,Chrom16_block2 ,Chrom16_block3 ,Chrom16_block4 ,Chrom16_block5 ,Chrom17_block1 ,Chrom18_block1 ,Chrom18_block2 ,Chrom19_block1 ,Chrom20_block1 ,Chrom20_block2 ,Chrom21_block1 ,Chrom21_block2 ,Chrom21_block3 ,Chrom22_block1 ,Chrom22_block2 ,Chrom23_block1 ,Chrom23_block2 ,Chrom23_block3 ,Chrom23_block4 ,Chrom24_block1 ,Chrom24_block2 ,Chrom25_block1 ,Chrom26_block1 ,Chrom26_block2 ,Chrom27_block1 ,Chrom28_block1 ,Chrom28_block2 ,Chrom29_block1 ,Chrom30_block1 ,Chrom31_block1 ,Chrom32_block1 ,Chrom33_block1 ,Chrom33_block2 ,Chrom34_block1 ,Chrom34_block2)

#874417618 total length of syntenic blocks
#902676312 total length in non-syntenic regions
#1777093930 total length of genome in chromosomes
#0.4920492% of genome is syntenic

write.csv(Tsch_syntenic_genes, file = "Tsch_syntenic_genes.csv")

Tsch_uni_ohnos <- read.csv("Tsch_universal_ohnos.csv")
Tsch_syntenic_ohnologs_all <-Tsch_uni_ohnos$Otsch.x%in%Tsch_syntenic_genes$Otsch
Tsch_redip <- read.csv("Othsch_all_redip_gene_ID.csv")
Tsch_syntenic_rediploidization <- Tsch_redip$Otsch.y%in%Tsch_syntenic_genes$Otsch

Tsch_syntenic_stuff <- cbind(Tsch_protins,Tsch_synteny_labels)
Tsch_downsampled_genes <- Tsch_syntenic_stuff[sample(1:nrow(Tsch_syntenic_stuff),8406,
                                                     replace=FALSE),]
table(Tsch_downsampled_genes$Tsch_synteny_labels)["TRUE"]
Tsch_uni_ohnos_synteny <- matrix(c(4713,8406,3403,8406),
                                 nrow = 2,
                                 dimnames = list(c("Hits","Total"),
                                                 c("Hits","Total")))
fisher.test(Tsch_uni_ohnos_synteny, alternative = 'two.sided')
#odds ratio 1.38492, p-value < 2.2e-16

Tsch_syntenic_stuff <- cbind(Tsch_protins,Tsch_synteny_labels)
Tsch_downsampled_genes <- Tsch_syntenic_stuff[sample(1:nrow(Tsch_syntenic_stuff),8406,
                                                     replace=FALSE),]
table(Tsch_downsampled_genes$Tsch_synteny_labels)["TRUE"]


Tsch_synteny_redip_all <- matrix(c(3753,10157,4074,10157),
                                 nrow = 2,
                                 dimnames = list(c("Hits","Total"),
                                                 c("Hits","Total")))
fisher.test(Tsch_synteny_redip_all, alternative = 'two.sided')
#odds ratio 0.9212107 p-value = 0.002116
Otsch_karotype <- read.csv("Otsh_karyotype.csv")
Otsch_synteny <- read.csv("Tsch_Chromosome_synteny_graph.csv")
ideogram(karyotype = Otsch_karotype, overlaid = Otsch_synteny,colorset1 = c('#30326F','#424EA0','#609FD4','#B899B8','#C3D8EB','#142511','#213F19','#71A45F','#F1EFC2','#A5A277','#856234','#CA9C49','#CC2E18','#44130C','#DCDEDD','#B8D7DC','#EDD37E','#DF8F46','#9A7356'))


##Sockeye
Nerka__protein <- read.csv("Nerka_proteins.csv")
Nerka__Chrom1_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG1' & Start >= 8875133 & Start <= 26014173)
Nerka__Chrom2_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG2' & Start >= 8183344 & Start <= 23791530)
Nerka__Chrom2_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG2' & Start >= 49170003 & Start <= 49277673)
Nerka__Chrom2_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG2' & Start >= 3388263 & Start <= 10590790)
Nerka__Chrom3_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG3' & Start <= 44018881 & Start <= 47419724)
Nerka__Chrom3_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG3' & Start >= 20702054 & Start <= 3098204)
Nerka__Chrom3_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG3' & Start >= 36962228 & Start <= 40863657)
Nerka__Chrom4_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG4' & Start >= 16367637 & Start <= 23271297)
Nerka__Chrom4_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG4' & Start >= 4425971 & Start <= 6315716)
Nerka__Chrom4_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG4' & Start >= 15566615 & Start <= 31460496)
Nerka__Chrom4_block4 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG4' & Start >= 42323712 & Start <= 47054889)
Nerka__Chrom5_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG5' & Start >= 28050751 & Start <= 28159500)
Nerka__Chrom5_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG5' & Start >= 21192248 & Start <= 21293744)
Nerka__Chrom5_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG5' & Start >= 19146745 & Start <= 2875895)
Nerka__Chrom5_block4 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG5' & Start >= 27727498 & Start <= 34959400)
Nerka__Chrom6_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG6' & Start <= 30739108 & Start >= 52308229)
Nerka__Chrom6_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG6' & Start >= 3751811 & Start <= 24207389)
Nerka__Chrom7_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG7' & Start >= 9847340 & Start <= 18152964)
Nerka__Chrom7_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG7' & Start >= 28866636 & Start <= 28866636)
Nerka__Chrom8_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG8' & Start >= 975969 & Start <= 28566556)
Nerka__Chrom9_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG9' & Start >= 14716626 & Start <= 37944406)
Nerka__Chrom9b_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG9b' & Start >= 5654237 & Start <= 26050151)
Nerka__Chrom10_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG10' & Start >= 6080284 & Start <= 16836502)
Nerka__Chrom10_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG10' & Start >= 30780591 & Start <= 50738676)
Nerka__Chrom11_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG11' & Start >= 13083787 & Start <= 17803775)
Nerka__Chrom11_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG11' & Start >= 17749322 & Start <= 22317139)
Nerka__Chrom12_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG12' & Start >= 6332437 & Start <= 15449295)
Nerka__Chrom12_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG12' & Start >= 17106716 & Start <= 27833023)
Nerka__Chrom13_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG13' & Start >= 9584209 & Start <= 29147175)
Nerka__Chrom13_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG13' & Start >= 33567412 & Start <= 52868197)
Nerka__Chrom14_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG14' & Start >= 26545554 & Start <= 43159080)
Nerka__Chrom14_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG14' & Start >= 1865952 & Start <= 20532312)
Nerka__Chrom14_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG14' & Start >= 31125614 & Start <= 31267238)
Nerka__Chrom15_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG15' & Start >= 18878716 & Start <= 19017846)
Nerka__Chrom15_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG15' & Start >= 5193297 & Start <= 11470738)
Nerka__Chrom15_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG15' & Start >= 27711894 & Start <= 31214401)
Nerka__Chrom15_block4 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG15' & Start >= 46577336 & Start <= 48187139)
Nerka__Chrom16_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG16' & Start >= 2694794 & Start <= 2884157)
Nerka__Chrom16_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG16' & Start >= 7671364 & Start <= 29744005)
Nerka__Chrom17_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG17' & Start >= 4188942 & Start <= 25910150)
Nerka__Chrom18_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG18' & Start >= 1519059 & Start <= 5547474)
Nerka__Chrom18_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG18' & Start >= 14244017 & Start <= 41663715)
Nerka__Chrom19_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG19' & Start >= 15462610 & Start <= 26169638)
Nerka__Chrom19_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG19' & Start >= 6703057 & Start <= 12891391)
Nerka__Chrom20_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG20' & Start >= 19839428 & Start <= 25344233)
Nerka__Chrom20_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG20' & Start >= 4536894 & Start <= 9355623)
Nerka__Chrom20_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG20' & Start >= 31790239 & Start <= 31966211)
Nerka__Chrom20_block4 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG20' & Start >= 1802995 & Start <= 1909081)
Nerka__Chrom20_block5 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG20' & Start >= 33828945 & Start <= 66646538)
Nerka__Chrom21_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG21' & Start >= 17830532 & Start <= 23450708)
Nerka__Chrom21_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG21' & Start >= 3343421 & Start <= 7172737)
Nerka__Chrom21_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG21' & Start >= 24805412 & Start <= 24926361)
Nerka__Chrom22_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG22' & Start >= 12955276 & Start <= 26060974)
Nerka__Chrom22_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG22' & Start >= 28710263 & Start <= 29916844)
Nerka__Chrom22_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG22' & Start >= 46829276 & Start <= 63922398)
Nerka__Chrom23_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG23' & Start >= 6900212 & Start <= 20389305)
Nerka__Chrom23_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG23' & Start >= 32664802 & Start <= 40012700)
Nerka__Chrom24_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG24' & Start >= 33960213 & Start <= 52088720)
Nerka__Chrom24_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG24' & Start >= 6714535 & Start <= 24529359)
Nerka__Chrom25_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG25' & Start >= 17595389 & Start <= 42644416)
Nerka__Chrom26_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG26' & Start >= 8488903 & Start <= 8594358)
Nerka__Chrom26_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG26' & Start >= 3165981 & Start <= 21813149)
Nerka__Chrom27_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG27' & Start >= 60888262 & Start <= 61008815)
Nerka__Chrom27_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG27' & Start >= 2959444 & Start <= 3122611)
Nerka__Chrom27_block3 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG27' & Start >= 10598147 & Start <= 18724301)
Nerka__Chrom27_block4 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG27' & Start >= 30193814 & Start <= 32160101)
Nerka__Chrom27_block5 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG27' & Start >= 44166820 & Start <= 61554015)
Nerka__Chrom28_block1 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG28' & Start >= 48370363 & Start <= 58398805)
Nerka__Chrom28_block2 <- filter(Nerka__protein, Nerka__Chromosome == 'linkage group LG28' & Start >= 26287200 & Start <= 40108092)

#708715686 total length of syntenic blocks
#1469280629 total length of genome in chromosomes
#0.4823556 percent of genomic is syntenic
Nerka_syntenic_genes <- rbind(Nerka__Chrom1_block1,Nerka__Chrom2_block1,Nerka__Chrom2_block2,Nerka__Chrom2_block3,Nerka__Chrom3_block1,Nerka__Chrom3_block2,Nerka__Chrom3_block3,Nerka__Chrom4_block1,Nerka__Chrom4_block2,Nerka__Chrom4_block3,Nerka__Chrom4_block4,Nerka__Chrom5_block1,Nerka__Chrom5_block2,Nerka__Chrom5_block3,Nerka__Chrom5_block4,Nerka__Chrom6_block1,Nerka__Chrom6_block2,Nerka__Chrom7_block1,Nerka__Chrom7_block2,Nerka__Chrom8_block1,Nerka__Chrom9_block1,Nerka__Chrom9b_block1,Nerka__Chrom10_block1,Nerka__Chrom10_block2,Nerka__Chrom11_block1,Nerka__Chrom11_block2,Nerka__Chrom12_block1,Nerka__Chrom12_block2,Nerka__Chrom13_block1,Nerka__Chrom13_block2,Nerka__Chrom14_block1,Nerka__Chrom14_block2,Nerka__Chrom14_block3,Nerka__Chrom15_block1,Nerka__Chrom15_block2,Nerka__Chrom15_block3,Nerka__Chrom15_block4,Nerka__Chrom16_block1,Nerka__Chrom16_block2,Nerka__Chrom17_block1,Nerka__Chrom18_block1,Nerka__Chrom18_block2,Nerka__Chrom19_block1,Nerka__Chrom19_block2,Nerka__Chrom20_block1,Nerka__Chrom20_block2,Nerka__Chrom20_block3,Nerka__Chrom20_block4,Nerka__Chrom20_block5,Nerka__Chrom21_block1,Nerka__Chrom21_block2,Nerka__Chrom21_block3,Nerka__Chrom22_block1,Nerka__Chrom22_block2,Nerka__Chrom22_block3,Nerka__Chrom23_block1,Nerka__Chrom23_block2,Nerka__Chrom24_block1,Nerka__Chrom24_block2,Nerka__Chrom25_block1,Nerka__Chrom26_block1,Nerka__Chrom26_block2,Nerka__Chrom27_block1,Nerka__Chrom27_block2,Nerka__Chrom27_block3,Nerka__Chrom27_block4,Nerka__Chrom27_block5,Nerka__Chrom28_block1,Nerka__Chrom28_block2)
write.csv(Nerka_syntenic_genes, file = "Onerka_syntenic_genes.csv")

Nerka_uni_ohnos <- read.csv("Nerka_universal_ohnos.csv")
Nerka_syntenic_ohnologs_all <-Nerka_uni_ohnos$Onerka.x%in%Nerka_syntenic_genes$Protein.product
table(Nerka_syntenic_ohnologs_all)["TRUE"]
Nerka_synteny_labels <- Nerka__protein$Protein.product%in%Nerka_syntenic_genes$Protein.product
table(Nerka_synteny_labels)["TRUE"]
Nerka_syntenic_stuff <- cbind(Nerka__protein,Nerka_synteny_labels)
Nerka_downsampled_genes <- Nerka_syntenic_stuff[sample(1:nrow(Nerka_syntenic_stuff),8406,
                                                       replace=FALSE),]
table(Nerka_downsampled_genes$Nerka_synteny_labels)["TRUE"]
Nerka_uni_ohnos_synteny <- matrix(c(4442,8406,2929,8406),
                                  nrow = 2,
                                  dimnames = list(c("Hits","Total"),
                                                  c("Hits","Total")))
fisher.test(Nerka_uni_ohnos_synteny, alternative = 'two.sided')
#odds ratio 1.516509, p-value < 2.2e-16

Onerka_redip <- read.csv("Onerka_redip_gene_ids.csv")
Onerka_syntenic_rediploidized <- Onerka_redip$Onerka.y%in%Nerka_syntenic_genes$Protein.product
table(Onerka_syntenic_rediploidized)["TRUE"]
Nerka_synteny_labels <- Nerka__protein$Protein.product%in%Nerka_syntenic_genes$Protein.product
table(Nerka_synteny_labels)["TRUE"]
Nerka_syntenic_stuff <- cbind(Nerka__protein,Nerka_synteny_labels)
Nerka_downsampled_genes <- Nerka_syntenic_stuff[sample(1:nrow(Nerka_syntenic_stuff),9969,
                                                       replace=FALSE),]
table(Nerka_downsampled_genes$Nerka_synteny_labels)["TRUE"]
#3305 syntenic rediploidized
#6664 non-syntenic rediploidized
#9813
Nerka_redip_location <- matrix(c(3305,9813,3589,9813),
                               nrow = 2,
                               dimnames = list(c("Hits","Total"),
                                               c("Hits","Total")))
fisher.test(Nerka_redip_location, alternative = 'two.sided')
#Odds ratio 0.9208733 p-value = 0.00328

Onerka_karotype <- read.csv("Onerka_karyotype.csv")
Onerka_synteny <- read.csv("Nerka_synteny_graph_input.csv")

head(Onerka_synteny)
ideogram(karyotype = Onerka_karotype, overlaid = Onerka_synteny,colorset1 = c('#30326F','#424EA0','#609FD4','#B899B8','#C3D8EB','#142511','#213F19','#71A45F','#F1EFC2','#A5A277','#856234','#CA9C49','#CC2E18','#44130C','#DCDEDD','#B8D7DC','#EDD37E','#DF8F46','#9A7356'))
convertSVG("chromosome.svg", device = "png")

#Rainbow Trout
Mykiss_protein <- read.csv("Mykiss_proteins.csv")
Mykiss_Chrom1_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 1' & Start >= 8879842 & Start <= 9022194)
Mykiss_Chrom1_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 1' & Start >= 30161228 & Start <= 32159116)
Mykiss_Chrom1_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 1' & Start >= 64763940 & Start <= 77560666)
Mykiss_Chrom1_block4 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 1' & Start >= 42125765 & Start <= 45357752)
Mykiss_Chrom2_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 2' & Start <= 23007051 & Start >= 23345653)
Mykiss_Chrom2_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 2' & Start >= 26026752 & Start <= 26214721)
Mykiss_Chrom2_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 2' & Start >= 34881661 & Start <= 36662480)
Mykiss_Chrom2_block4 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 2' & Start >= 53902359 & Start <= 85930792)
Mykiss_Chrom2_block5 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 2' & Start >= 37293193 & Start <= 42575378)
Mykiss_Chrom3_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 3' & Start >= 11578109 & Start <= 60706878)
Mykiss_Chrom4_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 4' & Start >= 24509335 & Start <= 24691893)
Mykiss_Chrom4_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 4' & Start >= 13590517 & Start <= 19192692)
Mykiss_Chrom5_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 5' & Start >= 33062937 & Start <= 43588482)
Mykiss_Chrom5_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 5' & Start >= 52403793 & Start <= 77340002)
Mykiss_Chrom5_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 5' & Start >= 13515068 & Start <= 29238350)
Mykiss_Chrom6_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 6' & Start <= 55403921 & Start >= 72659093)
Mykiss_Chrom6_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 6' & Start >= 18076575 & Start <= 39570937)
Mykiss_Chrom7_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 7' & Start >= 17025425 & Start <= 25755918)
Mykiss_Chrom7_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 7' & Start >= 37112312 & Start <= 54460731)
Mykiss_Chrom8_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 8' & Start >= 54412185 & Start <= 74176391)
Mykiss_Chrom8_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 8' & Start >= 23732904 & Start <= 47927715)
Mykiss_Chrom9_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 9' & Start >= 37566196 & Start <= 64794361)
Mykiss_Chrom9_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 9' & Start >= 23536372 & Start <= 26481903)
Mykiss_Chrom9_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 9' & Start >= 20641447 & Start <= 27471748)
Mykiss_Chrom10_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 10' & Start >= 27672654 & Start <= 66182670)
Mykiss_Chrom11_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 11' & Start >= 12099380 & Start <= 15976186)
Mykiss_Chrom11_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 11' & Start >= 34922495 & Start <= 38023790)
Mykiss_Chrom11_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 11' & Start >= 41686843 & Start <= 59403654)
Mykiss_Chrom11_block4 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 11' & Start >= 22213931 & Start <= 28034029)
Mykiss_Chrom12_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 12' & Start >= 15636624 & Start <= 54639233)
Mykiss_Chrom12_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 12' & Start <= 63558809 & Start >= 80446269)
Mykiss_Chrom13_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 13' & Start >= 65963868 & Start <= 66069927)
Mykiss_Chrom13_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 13' & Start >= 41987455 & Start <= 50403284)
Mykiss_Chrom14_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 14' & Start >= 19529551 & Start <= 42035889)
Mykiss_Chrom15_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 15' & Start >= 48609757 & Start <= 60293843)
Mykiss_Chrom15_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 15' & Start >= 20231954 & Start <= 31784729)
Mykiss_Chrom16_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 16' & Start >= 14881539 & Start <= 23977458)
Mykiss_Chrom16_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 16' & Start >= 40897561 & Start <= 48011115)
Mykiss_Chrom16_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 16' & Start >= 56857804 & Start <= 57819971)
Mykiss_Chrom16_block4 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 16' & Start >= 34355886 & Start <= 34532636)
Mykiss_Chrom16_block5 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 16' & Start >= 10195393 & Start <= 10332611)
Mykiss_Chrom17_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 17' & Start >= 25663667 & Start <= 77848756)
Mykiss_Chrom18_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 18' & Start >= 16013340 & Start <= 29521917)
Mykiss_Chrom18_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 18' & Start >= 32709615 & Start <= 51539263)
Mykiss_Chrom19_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 19' & Start >= 12154811 & Start <= 16817671)
Mykiss_Chrom19_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 19' & Start >= 25215660 & Start <= 58348522)
Mykiss_Chrom20_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 20' & Start >= 19093881 & Start <= 35420050)
Mykiss_Chrom20_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 20' & Start >= 11347875 & Start <= 16570402)
Mykiss_Chrom21_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 21' & Start >= 14414084 & Start <= 14520599)
Mykiss_Chrom21_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 21' & Start >= 15917750 & Start <= 16087667)
Mykiss_Chrom21_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 21' & Start >= 32202712 & Start <= 38745208)
Mykiss_Chrom22_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 22' & Start >= 32694843 & Start <= 32700451)
Mykiss_Chrom22_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 22' & Start >= 9673817 & Start <= 19695191)
Mykiss_Chrom22_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 22' & Start >= 24076693 & Start <= 24198516)
Mykiss_Chrom23_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 23' & Start >= 14455431 & Start <= 33975558)
Mykiss_Chrom24_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 24' & Start >= 1913755 & Start <= 27333237)
Mykiss_Chrom25_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 25' & Start >= 15071066 & Start <= 24492707)
Mykiss_Chrom26_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 26' & Start >= 14903564 & Start <= 35051725)
Mykiss_Chrom26_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 26' & Start >= 13120540 & Start <= 13221501)
Mykiss_Chrom27_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 27' & Start >= 32893838 & Start <= 33049762)
Mykiss_Chrom27_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 27' & Start >= 8623387 & Start <= 21463448)
Mykiss_Chrom28_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 28' & Start >= 12114412 & Start <= 28686192)
Mykiss_Chrom30_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 30' & Start >= 11836645 & Start <= 44774836)
Mykiss_Chrom31_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 31' & Start >= 17244472 & Start <= 17355043)
Mykiss_Chrom31_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 31' & Start >= 9542034 & Start <= 9674224)
Mykiss_Chrom31_block3 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 31' & Start >= 18596913 & Start <= 18733910)
Mykiss_Chrom31_block4 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 31' & Start >= 28060998 & Start <= 34958120)
Mykiss_Chrom32_block1 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 32' & Start >= 14414084 & Start <= 14520599)
Mykiss_Chrom32_block2 <- filter(Mykiss_protein, Mykiss_Chromosome == 'chromosome 32' & Start >= 16087667 & Start <= 38595107)

#851401738 total length of sytenic blocks
#2186890859 total length of genome in chromosomes
# 0.3893206 of nucleotides are syntenic
Mykiss_syntenic_genes <- rbind(Mykiss_Chrom1_block1,Mykiss_Chrom1_block2,Mykiss_Chrom1_block3,Mykiss_Chrom1_block4,Mykiss_Chrom2_block1,Mykiss_Chrom2_block2,Mykiss_Chrom2_block3,Mykiss_Chrom2_block4,Mykiss_Chrom2_block5,Mykiss_Chrom3_block1,Mykiss_Chrom4_block1,Mykiss_Chrom4_block2,Mykiss_Chrom5_block1,Mykiss_Chrom5_block2,Mykiss_Chrom5_block3,Mykiss_Chrom6_block1,Mykiss_Chrom6_block2,Mykiss_Chrom7_block1,Mykiss_Chrom7_block2,Mykiss_Chrom8_block1,Mykiss_Chrom8_block2,Mykiss_Chrom9_block1,Mykiss_Chrom9_block2,Mykiss_Chrom9_block3,Mykiss_Chrom10_block1,Mykiss_Chrom11_block1,Mykiss_Chrom11_block2,Mykiss_Chrom11_block3,Mykiss_Chrom11_block4,Mykiss_Chrom12_block1,Mykiss_Chrom12_block2,Mykiss_Chrom13_block1,Mykiss_Chrom13_block2,Mykiss_Chrom14_block1,Mykiss_Chrom15_block1,Mykiss_Chrom15_block2,Mykiss_Chrom16_block1,Mykiss_Chrom16_block2,Mykiss_Chrom16_block3,Mykiss_Chrom16_block4,Mykiss_Chrom16_block5,Mykiss_Chrom17_block1,Mykiss_Chrom18_block1,Mykiss_Chrom18_block2,Mykiss_Chrom19_block1,Mykiss_Chrom19_block2,Mykiss_Chrom20_block1,Mykiss_Chrom20_block2,Mykiss_Chrom21_block1,Mykiss_Chrom21_block2,Mykiss_Chrom21_block3,Mykiss_Chrom22_block1,Mykiss_Chrom22_block2,Mykiss_Chrom22_block3,Mykiss_Chrom23_block1,Mykiss_Chrom24_block2,Mykiss_Chrom25_block1,Mykiss_Chrom26_block1,Mykiss_Chrom26_block2,Mykiss_Chrom27_block1,Mykiss_Chrom27_block2,Mykiss_Chrom28_block1,Mykiss_Chrom30_block1,Mykiss_Chrom31_block1,Mykiss_Chrom31_block2,Mykiss_Chrom31_block3,Mykiss_Chrom31_block4,Mykiss_Chrom32_block1,Mykiss_Chrom32_block2)
write.csv(Mykiss_syntenic_genes, file = "Mykiss_syntenic_genes.csv")

Mykiss_uni_ohnos <- read.csv("Mykiss_universal_ohnos.csv")
Mykiss_syntenic_ohnologs_all <-Mykiss_uni_ohnos$Omykiss.x%in%Mykiss_syntenic_genes$Protein.product
table(Mykiss_syntenic_ohnologs_all)["TRUE"]
Mykiss_synteny_labels <- Mykiss_protein$Protein.product%in%Mykiss_syntenic_genes$Protein.product
Mykiss_syntenic_stuff <- cbind(Mykiss_protein,Mykiss_synteny_labels)
Mykiss_downsampled_genes <- Mykiss_syntenic_stuff[sample(1:nrow(Mykiss_syntenic_stuff),8406,
                                                         replace=FALSE),]
table(Mykiss_downsampled_genes$Mykiss_synteny_labels)["TRUE"]
Mykiss_uni_ohnos_synteny <- matrix(c(4801,8406,3414,8406),
                                   nrow = 2,
                                   dimnames = list(c("Hits","Total"),
                                                   c("Hits","Total")))
fisher.test(Mykiss_uni_ohnos_synteny, alternative = 'two.sided')
#odds ratio 1.406231, p-value < 2.2e-16


Mykiss_redip_genes <- read.csv("Omykiss_redip_gene_ids.csv")
#8990
Mykiss_syntenic_genes <- read.csv("Mykiss_syntenic_genes.csv")

Mykiss_syntenic_redips <- Mykiss_redip_genes$Omykiss.y%in%Mykiss_syntenic_genes$Protein.product
Mykiss_synteny_labels <- Mykiss_protein$Protein.product%in%Mykiss_syntenic_genes$Protein.product
Mykiss_syntenic_stuff <- cbind(Mykiss_protein,Mykiss_synteny_labels)
Mykiss_downsampled_genes <- Mykiss_syntenic_stuff[sample(1:nrow(Mykiss_syntenic_stuff),8990,
                                                         replace=FALSE),]
table(Mykiss_downsampled_genes$Mykiss_synteny_labels)["TRUE"]
#8605 repidploidized genes
#3708 syntenic redips
#5282 non-syntenic redips
Mykiss_redip_synteny <- matrix(c(3708,8990,3619,8990),
                               nrow = 2,
                               dimnames = list(c("Hits","Total"),
                                               c("Hits","Total")))
fisher.test(Mykiss_redip_synteny, alternative = 'two.sided')
#Odds ratio 1.024598 p-value = 0.38267

mykiss_karotype <- read.csv("Mykiss_karyotype.csv")
mykiss_synteny <- read.csv("Mykiss_synteny_graph_input.csv")
ideogram(karyotype = mykiss_karotype, overlaid = mykiss_synteny,colorset1 = c('#30326F','#424EA0','#609FD4','#B899B8','#C3D8EB','#142511','#213F19','#71A45F','#F1EFC2','#A5A277','#856234','#CA9C49','#CC2E18','#44130C','#DCDEDD','#B8D7DC','#EDD37E','#DF8F46','#9A7356'))




