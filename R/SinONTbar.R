library(B2B)
library(data.table)
library(ShortRead)
library(Biostrings)
library(DECIPHER)

#' Calculate two sequence similarity
#'
#' @param x please check parSeqSim input
#'
#' @return
#' @export
#'
#' @examples
SequenceSimilarity <- function (x)
{
  score <- rDNAse::parSeqSim(as.list(x), cores = 1, type = "local",
                             submat = "BLOSUM100")
  score <- reshape2::melt(score)$value
  if (mean(score) == 1) {
    1
  }
  else {
    mean(score[score != 1])
  }
}

#' extract barcode seq form ONT data according to 1 pattern
#'
#' @param read fastq read in DNAstring format
#' @param qual fastq quality
#' @param patta linker1
#' @param patta2 linker2
#' @param StartMismatch the mismatch base number threshold when searching for pattern in ONT reads using matchPattern
#' @param BarcodeLength Singleron  single barcode length
#' @param umilength Singleron  umi length
#' @param linkerlength Singleron single linker length
#' @param barcode_max_dis the maximum distance between two linkers, if there is no indels or delete, the distance is BarcodeLength 8.
#'
#' @return Barcode + umi sequence in DNAstringset format extracted from ONT data.
#' @export
#'
#' @examples
ONTBarcodeSeq <- function (read, qual, patta, patta2, StartMismatch = 5,
                        BarcodeLength = 8, umilength = 12, linkerlength = 16, barcode_max_dis = 12)
{
  # group1 patt1
    pattaRC <- Biostrings::reverseComplement(patta)
    mismatch <- StartMismatch
    mP <- Biostrings::matchPattern(pattern = patta, subject = read,
        max.mismatch = mismatch, with.indels = TRUE)
    mPRC <- Biostrings::matchPattern(pattern = pattaRC, subject = read,
        max.mismatch = mismatch, with.indels = TRUE)
  # group1 patt2
    pattaRC2 <- Biostrings::reverseComplement(patta2)
    mismatch <- StartMismatch
    mP2 <- Biostrings::matchPattern(pattern = patta2, subject = read,
        max.mismatch = mismatch, with.indels = TRUE)
    mPRC2 <- Biostrings::matchPattern(pattern = pattaRC2, subject = read,
        max.mismatch = mismatch, with.indels = TRUE)

  # judge strand
    if (length(mP) != 0 & length(mPRC) != 0) {

      mPSco <-  unlist(lapply(mP, function(x){sco <- SequenceSimilarity(as.character(
        c(Biostrings::DNAStringSet(patta), Biostrings::DNAStringSet(x))))}))
      mPLin <- mP[which.max(mPSco)]
      mPSco <-  max(mPSco)

      mPRCSco = unlist(lapply(mPRC, function(x){sco <- SequenceSimilarity(as.character(
        c(Biostrings::DNAStringSet(pattaRC), Biostrings::DNAStringSet(x))))}))
      mPRCLin =  mPRC[which.max(mPRCSco)]
      mPRCSco = max(mPRCSco)

      mPSco2 <-  unlist(lapply(mP2, function(x){sco <- SequenceSimilarity(as.character(
        c(Biostrings::DNAStringSet(patta2), Biostrings::DNAStringSet(x))))}))
      mPLin2 <- mP2[which.max(mPSco2)]
      mPSco2 <-  max(mPSco2)

      mPRCSco2 = unlist(lapply(mPRC2, function(x){sco <- SequenceSimilarity(as.character(
        c(Biostrings::DNAStringSet(pattaRC2), Biostrings::DNAStringSet(x))))}))
      mPRCLin2 =  mPRC2[which.max(mPRCSco2)]
      mPRCSco2 = max(mPRCSco2)

      mPSco = mean(mPSco, mPSco2)
      mPRCSco2 = mean(mPRCSco, mPRCSco2)

      dis = abs(start(mPLin2) - end(mPLin) - 1)
      disRS = abs(start(mPRCLin) - end(mPRCLin2) - 1)
      if(min(dis, disRS) < barcode_max_dis){
        if (mPSco == mPRCSco) {
          if(abs(dis - 8) < abs(disRS - 8)) {
            TR1 <- mP
            Orien <- 1
          }else{
            TR1 <- mPRC
            Orien <- 2
          }
        } else if (mPSco > mPRCSco) {
                TR1 <- mP
                Orien <- 1
            }
          else {
                TR1 <- mPRC
                Orien <- 2
            }
      } else {
      return(NA)}
    } else if (length(mP) != 0 & length(mPRC) == 0) {
        TR1 <- mP
        Orien <- 1
    } else if (length(mP) == 0 & length(mPRC) != 0) {
        TR1 <- mPRC
        Orien <- 2
    } else if (length(mP) == 0 & length(mPRC) == 0) {
        return(NA)
    }
  # extract seq
    if (Orien == 1) {
      start(TR1) <- start(TR1) - BarcodeLength
      end(TR1) <- end(TR1) + BarcodeLength + linkerlength + BarcodeLength + 1 + umilength
      TR1 <- TR1[IRanges::width(Biostrings::DNAStringSet(TR1)) ==
            IRanges::width(TR1)]
      if(length(TR1) > 1) {
        TR2 = TR1
        start(TR2)  = start(TR1) + BarcodeLength + linkerlength + BarcodeLength
        end(TR2) = end(TR2) - BarcodeLength - umilength
        TR2 <- TR2[IRanges::width(Biostrings::DNAStringSet(TR2)) ==
            IRanges::width(TR2)]
        sco = unlist(lapply(TR2, function(x){sco <- SequenceSimilarity(as.character(c(Biostrings::DNAStringSet(patta2),
                                                 Biostrings::DNAStringSet(x))))
                  return(sco)}))
        TR1 = TR1[which.max(sco)]
        Barcode <- Biostrings::DNAStringSet(TR1)
        names(Barcode) <- 1
      }else if (length(TR1) == 1){
        TR1 <- TR1
        Barcode <- Biostrings::DNAStringSet(TR1)
        names(Barcode) <- 1
      } else if (length(TR1) == 0){
        return(NA)
      }
    }else if(Orien == 2){
      start(TR1) <- start(TR1) - BarcodeLength - linkerlength - BarcodeLength - 1 - umilength
      end(TR1) <- end(TR1) + BarcodeLength
      TR1 <- TR1[IRanges::width(Biostrings::DNAStringSet(TR1)) ==
            IRanges::width(TR1)]
      if(length(TR1) > 1) {
        TR2 = TR1
        start(TR2)  = start(TR1) + umilength + 1 + BarcodeLength
        end(TR2) = end(TR2) -  BarcodeLength - linkerlength - BarcodeLength
        TR2 <- TR2[IRanges::width(Biostrings::DNAStringSet(TR2)) ==
            IRanges::width(TR2)]
        pattaRC2 <- Biostrings::reverseComplement(patta2)
        sco = unlist(lapply(TR2, function(x){sco <- SequenceSimilarity(as.character(c(Biostrings::DNAStringSet(pattaRC2),
                                                 Biostrings::DNAStringSet(x))))
                  return(sco)}))
        TR1 = TR1[which.max(sco)]
        Barcode <- Biostrings::DNAStringSet(TR1)
        Barcode <- Biostrings::reverseComplement(Barcode)
        names(Barcode) <- 2
      } else if (length(TR1) == 1){
        TR1 <- TR1
        Barcode <- Biostrings::DNAStringSet(TR1)
        Barcode <- Biostrings::reverseComplement(Barcode)
        names(Barcode) <- 2
      } else if (length(TR1) == 0){
        return(NA)
      }
    }
    return(Barcode)
}



#' extract barcode seq form ONT data according to 4 pattern
#'
#' @param ONTfastq ONTfastq site, for example: /mnt/data1/zhangdan/Proj/Mouse_iso_atlas/scnapbar/sc_NGScDNA_C57_F_2m_c-kit+/data/nanopore_small.fq
#'
#' @return ONT barcode sequence
#' @export
#'
#' @examples
PattafourSeq <- function(ONTfastq){
  # linker pattern
  grp1_patta <- DNAString("TCGGTGACAGCCATAT")
  grp1_patta2 <- DNAString("CGTAGTCAGAAGCTGA")

  grp2_patta <- DNAString("ATCCACGTGCTTGAGA")
  grp2_patta2 <- DNAString("TCAGCATGCGGCTACG")

  grp3_patta <- DNAString("CGAACATGTAGGTCTC")
  grp3_patta2 <- DNAString("GACTACGTATTAGCAT")

  grp4_patta <- DNAString("GATTGTCACTAACGCG")
  grp4_patta2 <- DNAString("ATGCTGACTCCTAGTC")

  pattn = list(grp1_patta,grp1_patta2, grp2_patta, grp2_patta2,grp3_patta,grp3_patta2, grp4_patta, grp4_patta2)

  grp1_pattn_c = DNAString(paste0("NNNNNNNN", pattn[[1]], "NNNNNNNN", pattn[[2]], "NNNNNNNNC", "NNNNNNNNNNNN"))
  grp2_pattn_c = DNAString(paste0("NNNNNNNN", pattn[[3]], "NNNNNNNN", pattn[[4]], "NNNNNNNNC", "NNNNNNNNNNNN"))
  grp3_pattn_c = DNAString(paste0("NNNNNNNN", pattn[[5]], "NNNNNNNN", pattn[[6]], "NNNNNNNNC", "NNNNNNNNNNNN"))
  grp4_pattn_c = DNAString(paste0("NNNNNNNN", pattn[[7]], "NNNNNNNN", pattn[[8]], "NNNNNNNNC", "NNNNNNNNNNNN"))

  pattn_c = list(grp1_pattn_c, grp2_pattn_c, grp3_pattn_c, grp4_pattn_c)

  # ONTfastq precess
  fastq <- readFastq(ONTfastq)
  fqqual <- quality(fastq)
  IDTab <- do.call(rbind, base::strsplit(as.character(ShortRead::id(fastq)), "[ =\\|]"))
  IDTab <- as.data.table(IDTab)
  fastq <- sread(fastq)
  names(fastq) <- IDTab$V1

  # run ONTBarcodeSeq for 4 patterns
  mclapply(seq_along(fastq), function(x) {
    Bs_tmp =  lapply(c(1,3,5,7), function(i){
      patta1 = pattn[[i]]
      patta2 =  pattn[[i + 1]]
      Bs_tmp =  ONTBarcodeSeq(read = fastq[[x]], qual = fqqual[[x]], patta = patta1, patta2 = patta2,
                           StartMismatch = 5, BarcodeLength = 8, umilength = 12,
                           linkerlength = 16, barcode_max_dis = 12)
      return(Bs_tmp)
    })
    length =  unlist(lapply(Bs_tmp,function(x){nchar(as.character(x))}))
    #sco = unlist(lapply(1:2, function(x){sco <- SequenceSimilarity(as.character(c(Biostrings::DNAStringSet(pattn_c[[x]]),Biostrings::DNAStringSet(Bs_tmp[[x]]))))}))
    if (any(!is.na(length))){
      sco = unlist(lapply(1:4, function(x){sco <- adist(as.character(pattn_c[[x]]),as.character(Bs_tmp[[x]]), partial = TRUE)}))
      Bs = Bs_tmp[[which.min(sco)]]
      print(paste0("now read number is ", x))
      return(Bs)
    }
  }, mc.cores = 1) -> Bs

  BsTab <- data.table(Read = IDTab$V1, Seq = as.character((lapply(Bs, function(x){as.character(x)}))), Strand = as.character(lapply(Bs, function(x){names(x)})))
  BsTab$Strand <- gsub("NULL", NA, BsTab$Strand)
  BsTab <- as.data.frame(BsTab[!is.na(BsTab$Strand), ])
  BCs <-DNAStringSet(BsTab$Seq)
  names(BCs) <- BsTab$Read
  BCs <- BCs[width(BCs) >= 65]
  BCs
  return(BCs)
}



#' Match the ONT barcodes with Illunima NGS barcodes
#'
#' @param ONTfastq ONTfastq site, eg. /mnt/data1/zhangdan/Proj/Mouse_iso_atlas/scnapbar/sc_NGScDNA_C57_F_2m_c-kit+/data/nanopore_small.fq
#' @param Sinbarcode Singleron barcode site, eg. /mnt/data1/zhangdan/Proj/Mouse_iso_atlas/scnapbar/sc_NGScDNA_C57_F_2m_c-kit+/data/barcodes.tsv.gz
#' @param MaxMisMatchvalue barcode assignment MaxMisMatch, we got 57 base in total barcode, and 10 MaxMisMatch is set as default.
#'
#' @return Singleron NGS barcode and ONT read id assignment table
#' @export
#'
#' @examples
BarcodeAssign <- function(ONTfastq, Sinbarcode, MaxMisMatchvalue = 10){
  BCs = PattafourSeq(ONTfastq = ONTfastq)
  NGS <- fread(Sinbarcode, header = FALSE)
  ONTSB <- subseq(BCs, 1, 57)
  NGSSB <- DNAStringSet(unique(NGS$V1))
  MaxMisMatchvalue = MaxMisMatchvalue
  ONT2NGS <- mclapply(seq_along(ONTSB), function(i) B2B:::MatchSB2(ONT = ONTSB[[i]], SB = NGSSB, MaxMisMatch = MaxMisMatchvalue), mc.cores = 1)
  ONT2NGS <- as.data.frame(data.table(Read = names(ONTSB), SB = mapply(as.character, ONT2NGS)))
  ONT2NGS <- ONT2NGS[!is.na(ONT2NGS$SB), ]

  R2B0 <- merge(BsTab, ONT2NGS, by = "Read")
  return(R2B0)
}






