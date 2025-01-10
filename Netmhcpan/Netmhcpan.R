### fasta-HLA pseudosequence-format
    #mhci info come from netmhcpan4.1 training data; mhcii info come from netmhciipan4.3 training data
    load(paste0(path_test,"TEST.Rdata"))
    hlamhci01<-"/data1/wuguojia/data/mhc_benchmark/database/HLAallele_sequence/MHC_pseudo.dat"
    hlamhci02<-"/data1/wuguojia/data/mhc_benchmark/database/HLAallele_sequence/MHC_pseudo.dat"
    hlamhcii01<-"/data1/wuguojia/data/mhc_benchmark/database/HLAallele_sequence/pseudosequence.2023.dat"
    hlamhcii02<-"/data1/wuguojia/data/mhc_benchmark/database/HLAallele_sequence/pseudosequence.2016.all.X.dat"
    seq1<-rbind(read.table(hlamhci01, header = FALSE, sep = "", stringsAsFactors = FALSE),read.table(hlamhci02, header = FALSE, sep = "", stringsAsFactors = FALSE))
    seq2<-rbind(read.table(hlamhcii01, header = FALSE, sep = "", stringsAsFactors = FALSE),read.table(hlamhcii02, header = FALSE, sep = "", stringsAsFactors = FALSE))
    #mhci match sequences
    seq1$V1<-gsub("[^a-zA-Z0-9]", "", seq1$V1)
    seq1<-seq1 %>% distinct() %>% filter(grepl("^HLA",V1)) %>% distinct(V1,.keep_all=TRUE)
    setdiff(gsub("[^a-zA-Z0-9]","",unique(test_mhci$hlatype)),seq1$V1)
    #mhcii match sequences
    seq2$V1<-gsub("[^a-zA-Z0-9]", "", seq2$V1)
    seq2$V1 <- ifelse(startsWith(seq2$V1, "D"),paste0("HLA", seq2$V1), seq2$V1)
    seq2<-seq2 %>% distinct() %>% distinct(V1,.keep_all=TRUE)
    setdiff(gsub("[^a-zA-Z0-9]","",unique(test_mhcii$hlatype)),seq2$V1)
    #write fasta files
    hlaseq<-rbind(seq1,seq2)
    test_mhci$V1<-gsub("[^a-zA-Z0-9]", "", test_mhci$hlatype)
    test_mhcii$V1<-gsub("[^a-zA-Z0-9]", "", test_mhcii$hlatype)
    test_mhci<-left_join(test_mhci,hlaseq,by="V1")
    test_mhcii<-left_join(test_mhcii,hlaseq,by="V1")
    for(i in c("_mhci","_mhcii")){
        temp_tes <- get(paste0("test",i))
        temp_sel <- get(paste0("h2ltest",i))
        for(m in 1:nrow(temp_sel)){
            for(n in 1:(ncol(temp_sel)-1)){
                hla <- rownames(temp_sel)[m]
                len <- colnames(temp_sel)[n]
                num <- temp_sel[hla,len]
                if(num != 0){
                    sub <- temp_tes %>% filter(hlatype == hla & antigen_peptide_length == as.numeric(len))
                    hla <- gsub("/","-",hla)#linux don't allow '/' 
                    f <- file(paste0(path_HLA_fasta,hla,"_",len,".fasta"), open = "w")
                    for (k in 1:nrow(sub)) {
                        description <- paste0(">", sub$hlatype[k], "\n")
                        cat(description, file = f)
                        cat(sub$V2[k], file = f)
                        if (k < nrow(sub)) {cat("\n", file = f)}
                    }
                    close(f)
                }
            }
        }
    }
    files <- list.files("/data1/wuguojia/data/mhc_benchmark/testbase/HLA_fasta", recursive = TRUE, full.names = TRUE,pattern = "\\.fasta$")
    archive_mhci <- files[!grepl("HLA-D",files)]
    archive_mhcii <- files[grepl("HLA-D",files)]
    write.table(archive_mhci,file=paste0(path_HLA_fasta,"archive_mhci.txt"),col.names = F,row.names = F,quote = F)
    write.table(archive_mhcii,file=paste0(path_HLA_fasta,"archive_mhcii.txt"),col.names = F,row.names = F,quote = F)