## analyze TCGA junctions from TCGA based sequencing data
## TCGA Firehose pipeline

OutSplice_TCGA<-function(junction, gene_expr, rawcounts, output_file_prefix, dir, filterSex=T, genome = 'Homo.sapiens', annotation = 'org.Hs.eg.db', TxDb = 'TxDb.Hsapiens.UCSC.hg19.knownGene', offsets_value = 0.00001, correction_setting='fdr', p_value=0.05){

  date<-Sys.Date()

  ## this is non-log transformed data and includes pheno (RAW, in RPM)
  print('Loading data')
  all.junc<-read.table(file=junction, sep='\t', header=T, stringsAsFactors = F)
  all.gene_expr<-read.table(file=gene_expr, header=T, row.names=1, sep="\t")
  rawcounts<-read.table(file=rawcounts, sep='\t', header=T, row.names=1, stringsAsFactors = F)

  ##load info
  #source(file=paste0(dir,"OGSAfunctionwFisher.R"))

  ## Load needed libraries
  suppressPackageStartupMessages(library(genome, character.only = TRUE))
  suppressPackageStartupMessages(library(annotation, character.only = TRUE))
  suppressPackageStartupMessages(library('GenomicRanges'))
  suppressPackageStartupMessages(library('limma'))
  suppressPackageStartupMessages(library('Repitools'))
  suppressPackageStartupMessages(library(TxDb, character.only = TRUE))
  colnames(rawcounts)<-substr(colnames(rawcounts), 1,15)
  colnames(rawcounts)<-gsub('\\.', "-", colnames(rawcounts))

  #remove header; rename row names and junctions
  all.junc<-all.junc[-1,] # remove header
  all.gene_expr<-all.gene_expr[-1,] # remove header

  all.samples<-colnames(all.junc)
  all.samples<-substr(all.samples, 1, 15)
  all.samples<-gsub('\\.', '-', all.samples)
  junction.names<-all.junc[,1]
  #remove duplicates
  all.junc<-all.junc[!duplicated(junction.names),]
  junction.names<-all.junc[,1]
  j<-paste0(sapply(strsplit(junction.names,split=":"), function(x){x[[1]]}), ":", sapply(strsplit(junction.names,split=":"), function(x){x[[2]]}), "-",sapply(strsplit(junction.names,split=":"), function(x){x[[4]]}))

  rownames(all.junc)<-j
  colnames(all.junc)<-all.samples
  #remove first column
  all.junc<-all.junc[,-1]

  expression.samples<-substr(colnames(all.gene_expr), 1, 15)
  expression.samples<-gsub('\\.', '-', expression.samples)
  colnames(all.gene_expr)<-expression.samples

  #infer phenotype from sample names
  all.samples<-intersect(colnames(all.junc), colnames(all.gene_expr))
  all.samples<-intersect(all.samples, colnames(rawcounts))
  pheno<-substr(all.samples, 14,15)
  names(pheno)<-all.samples
  # include only primary tumor 01 or normal 11
  pheno<-pheno[pheno=='01'|pheno=="11"]
  pheno[pheno=='01']<-"Tumor"
  pheno[pheno=='11']<-"Normal"
  print(table(pheno))

  if (sum(pheno=="Normal")<10){
    stop('Too few normal samples')
  }
  #subset only primary tumor 01 or normal 11
  all.samples<-names(pheno)
  all.junc<-all.junc[,all.samples]
  all.gene_expr<-all.gene_expr[,all.samples]
  rawcounts<-rawcounts[,all.samples]
  pheno<-pheno[all.samples]

  ##change from char to numeric
  n<-sapply(all.gene_expr, as.numeric)
  rownames(n)<-rownames(all.gene_expr)
  all.gene_expr<-n
  n<-sapply(all.junc, as.numeric)
  rownames(n)<-rownames(all.junc)
  all.junc<-n
  remove(n)

  ## get data in RPM
  print("convert to RPM")
  totalrawcount<-colSums(rawcounts)
  names(totalrawcount)<-colnames(rawcounts)
  junc.RPM<-apply(all.junc[,all.samples], 1, function(x){x/totalrawcount[all.samples]*1000000})
  junc.RPM<-t(junc.RPM)
  #print(all.junc[1:10, 1:10]/junc.RPM[1:10, 1:10])

  print("filter the putative junctions")
  ### filter genes by cut off of
  fcCutoff=10
  Cutoff.ratio<-1-(1/fcCutoff)
  # enforce overall Fold change
  logFC <- apply(junc.RPM,1,function(x){(max(x)-min(x))/max(x)}) > Cutoff.ratio
  junc.RPM <- junc.RPM[which(logFC),]

  rm(logFC)
  gc()
  ###### NEW FILTER: See if tumors have any outliers ###########################
  ##PHENO should have 'Normal' or 'Tumor' calls where Tumor ==1, Normal ==0, and names of each sample associated
  PHENO<-pheno=='Tumor'
  PHENO<-as.numeric(PHENO)
  names(PHENO)<-names(pheno)

  print("run the ogsa function for pre filtering")
  ## get function
  test2<-dotheogsa(Sample.data=junc.RPM, PHENO=PHENO, offsets=0.1, dir = dir)
  has.outliers<-test2[,"outRankTumor1"]>1|test2[,"outRankTumor2"]>1
  junc.RPM<-junc.RPM[has.outliers,]

  ### Filter all genes on the X and Y chromosomes

  if (filterSex) {
    junc.RPM <- junc.RPM[grep('chr[XY]',row.names(junc.RPM),value=T,invert=T),]
  }
  ##############################################################
  print("get the genomic information for all the junctions")

  # create GenomicRanges object for junctions
  chr <- sapply(strsplit(row.names(junc.RPM),split=":"), function(x){x[[1]]})
  start <- as.numeric(sapply(strsplit(row.names(junc.RPM),split="[:-]"),
                             function(x){x[[2]]}))
  end <- as.numeric(sapply(strsplit(row.names(junc.RPM),split="[:-]"),
                           function(x){x[[3]]}))


  geneAnnot <- GRanges(seqnames=Rle(chr),
                       IRanges(start=start, end=end))
  names(geneAnnot) <- row.names(junc.RPM)

  geneAnnotAll <- GRanges(seqnames=Rle(chr),
                          IRanges(start=start, end=end))
  names(geneAnnotAll) <- row.names(junc.RPM)

  # get GenomicRanges object with genes for whole genome
  gn <- genes(get(TxDb))

  gSymbol <- select(get(annotation),keys=as.character(gn$gene_id),
                    columns=c('SYMBOL'),keytype='ENTREZID')
  gn$SYMBOL <- gSymbol$SYMBOL
  gn$ENTREZID <- gSymbol$ENTREZID
  rm(gSymbol, chr, end, start)
  gc()

  # find symbols and ENTREZID for our junctions
  overlap <- findOverlaps(geneAnnot,gn)

  geneSYMBOLS <- tapply(gn$SYMBOL[subjectHits(overlap)],
                        queryHits(overlap),paste,collapse=';')

  # add to genome ranges object
  geneAnnot$SYMBOL <- NA
  geneAnnot$SYMBOL[as.numeric(names(geneSYMBOLS))] <- geneSYMBOLS

  geneENTREZID <- tapply(gn$ENTREZID[subjectHits(overlap)],
                         queryHits(overlap),paste,collapse=';')

  geneAnnot$ENTREZID <- NA
  geneAnnot$ENTREZID[as.numeric(names(geneENTREZID))] <- geneENTREZID

  # get known exons
  en <- exons(get(TxDb))

  print("filtering for putative junctions, skipping, insertion and deletion events based on known exons.")
  # annotate skipping events as those which overlap with more than two exons
  exo <- findOverlaps(GRanges(seqnames(geneAnnot),
                              ranges=IRanges(start=start(geneAnnot)+1,
                                             end=end(geneAnnot)-1),
                              strand=strand(geneAnnot)),en,
                      type='within')
  exoCount <- table(factor(queryHits(exo))) > 0
  geneAnnot$skipping <- F
  geneAnnot$skipping[as.numeric(names(exoCount))] <- exoCount

  # insertion events start or end outside of known exons
  geneAnnot.start <- GRanges(seqnames(geneAnnot),
                             ranges=IRanges(start=start(geneAnnot),
                                            end=start(geneAnnot)),
                             strand=strand(geneAnnot))
  geneAnnot.end <- GRanges(seqnames(geneAnnot),
                           ranges=IRanges(start=end(geneAnnot),
                                          end=end(geneAnnot)),
                           strand=strand(geneAnnot))

  geneAnnot$insertions <- !overlapsAny(geneAnnot.start, en) |
    !overlapsAny(geneAnnot.end, en)

  # deletion events occur w/in an exon but not at its start or end
  en.start <- GRanges(seqnames(en),
                      ranges=IRanges(start=start(en),
                                     end=start(en)),
                      strand=strand(en))
  en.end <- GRanges(seqnames(en),
                    ranges=IRanges(start=end(en),
                                   end=end(en)),
                    strand=strand(en))
  geneAnnot$deletions <-
    (!overlapsAny(geneAnnot.end,en.start) & overlapsAny(geneAnnot.end,en)) |
    (!overlapsAny(geneAnnot.start,en.end) & overlapsAny(geneAnnot.start,en))

  filterToEvents<-T

  print("deletions")
  sum(geneAnnot$deletions)
  print("insertions")
  sum(geneAnnot$insertions)
  print("skipping")
  sum(geneAnnot$skipping)


  if (filterToEvents) {
    geneAnnot <- geneAnnot[apply(cbind(geneAnnot$deletions,
                                       geneAnnot$insertions,
                                       geneAnnot$skipping),1,any),]
    junc.RPM <- junc.RPM[names(geneAnnot),]
  }


  ##################################################################
  print("remove all that map to 'NA' no gene name, and assign gene expression from gene_expr")
  junc.RPM<-junc.RPM[!is.na(geneAnnot$SYMBOL),]
  geneAnnot<-geneAnnot[row.names(junc.RPM),]

  ## remove row names for unknown genes containing "?" unknown genes

  all.gene_expr <- all.gene_expr[!grepl('\\?', row.names(all.gene_expr)),]

  ## get gene names
  gene_exprgenenames<-sapply(strsplit(row.names(all.gene_expr), split ='\\|'), function(x){x[1]})
  gene_exprEntrezID<-sapply(strsplit(row.names(all.gene_expr), split ='\\|'), function(x){x[2]})

  print("align with gene_expr data")
  # ## get gene names


  ## collect gene_expr values
  # initialize matrix of values for each junction
  junctionGenegene_expr <- matrix(0, nrow=length(geneAnnot), ncol=ncol(all.gene_expr), dimnames=list(names(geneAnnot), colnames(all.gene_expr)))

  geneAnnot$ENTREZID -> genes2Junc_ENTREZ
  names(genes2Junc_ENTREZ)<-names(geneAnnot)
  print("shows how many junctions aligned to a single gene")
  length(grep(genes2Junc_ENTREZ,pattern=';',invert=T,value=T))
  ## this selects just the first gene that each junction aligns to
  genes2Junc_ENTREZ<-sapply(strsplit(genes2Junc_ENTREZ, ";"), function(x){x[1]})

  ## fill in the matrix
  no.gene_expr<-vector(length=length(genes2Junc_ENTREZ))
  for (g in 1:length(genes2Junc_ENTREZ)){
    if (!isTRUE(intersect(gene_exprEntrezID, genes2Junc_ENTREZ[g])>0)){
      ## for gene where there is no gene_expr data, skip it
      no.gene_expr[g]<-T
      next
    }
    junctionGenegene_expr[names(genes2Junc_ENTREZ)[g],]<-all.gene_expr[which(gene_exprEntrezID==genes2Junc_ENTREZ[g]),]

  }

  ################################################################
  ## Remove genes without gene_expr data, not in true genes
  junctionGenegene_expr[!no.gene_expr,]->junctionGenegene_expr
  junc.RPM<-junc.RPM[row.names(junctionGenegene_expr),]


  print("subset removing any genes without normalization")
  junc.RPM2<-junc.RPM[!(apply(junctionGenegene_expr, 1,sum))==0,]
  dim(junc.RPM2)
  # ## filters junctions
  junc.RPM.original<-junc.RPM
  junc.RPM<-junc.RPM2



  print("Perform normalization using gene_expr values")

  ## Turn all the zeros in junctionGenegene_expr into 1.
  ## replace ALL zeros with 1.  Then when you divide it does not do anything.
  junctionGenegene_expr<-junctionGenegene_expr[row.names(junc.RPM2),]
  junctionGenegene_expr2<-junctionGenegene_expr
  junctionGenegene_expr2[junctionGenegene_expr==0]<-1
  junc.RPM.norm<-junc.RPM2
  junc.RPM.norm<-junc.RPM2/junctionGenegene_expr2

  ### perform outlier analysis with OGSA
  ############### Use OGSA for outlier ranking #################

  ##PHENO should have 'Normal' or 'Tumor' calls where Tumor ==1, Normal ==0, and names of each sample associated
  PHENO<-pheno=='Tumor'
  PHENO<-as.numeric(PHENO)
  names(PHENO)<-names(pheno)

  print("run the ogsa function")

  FisherAnalyses=dotheogsa(Sample.data=junc.RPM.norm, PHENO=PHENO, offsets=offsets_value, Fisher=T, correction='fdr', dir = dir)
  ## use default offset=0.001 for normalized data
  # create lists of top genes
  topgenelist10.ogsa=FisherAnalyses[,'var1']

  topgenelist90.ogsa=FisherAnalyses[,'var2']


  ##Use Fisher to subset candidates that have Fisher test p value <0.05, or cutoff can be adjusted
  ## Underexpression in tumors
  toplist10<-FisherAnalyses[topgenelist10.ogsa,]
  if (sum(toplist10[,"FisherP1"]<p_value)>0){
    toplist10<-toplist10[toplist10[,"FisherP1"]<p_value,]
    toplist10<-toplist10[,"FisherP1"]

  }else{
    toplist10<-vector()
    print('no under expression outliers')
  }
  ## Overexpression in tumors
  toplist90<-FisherAnalyses[topgenelist90.ogsa,]
  if (sum(toplist90[,"FisherP2"]<p_value)>0){
    toplist90<-toplist90[toplist90[,"FisherP2"]<p_value,]
    toplist90<-toplist90[,"FisherP2"]

  }else{
    toplist90<-vector()
    print('no over expression outliers')
  }
  junctions<-c(names(toplist90), names(toplist10))
  pvalues<-c(toplist90, toplist10)

  ### get outier calls##
  junc.Outliers<-dotheogsa(Sample.data=junc.RPM.norm, PHENO=PHENO, offsets=offsets_value, Fisher=T, correction=correction_setting, outliers=T, dir = dir)

  ## Calculate median normalized value (of normals)
  ## for gene_expr
  #NORM.gene_expr<-junctionGenegene_expr[,pheno == 'Normal']
  ## get median of normal expression
  #NORM.gene_expr<-apply(NORM.gene_expr, 1, median)

  ## get junction expression in normal samples only
  #NORM.RPM<-junc.RPM[,pheno == 'Normal']
  ## get median of normal expression
  #NORM.RPM<-apply(NORM.RPM, 1, median)

  ## median of normal expression within the normalized data
  NORM.gene_expr.norm<-apply(junc.RPM.norm[,pheno == 'Normal'],1,median)

  ## aggregate the data
  junc.RPM<-junc.RPM[junctions,]
  gene_expr<-junctionGenegene_expr[junctions,]
  junc.RPM.norm<-junc.RPM.norm[junctions,]
  geneAnnot<-geneAnnot[junctions]
  ASE.type<-cbind(geneAnnot$skipping, geneAnnot$insertions, geneAnnot$deletions)
  colnames(ASE.type)<-c("skipping", "insertions", "deletions")
  row.names(ASE.type)<-junctions

  ## Calculate Splicing Burden
  #source(file=paste0(dir,"SpliceBurdenfunction.R"))
  splice_burden <- CalcBurden(junc.Outliers, FisherAnalyses, p_value)

  ## save output file
  save(junc.RPM, gene_expr, junc.RPM.norm, pvalues, pheno, FisherAnalyses, geneAnnot, ASE.type, NORM.gene_expr.norm, junc.Outliers, splice_burden, file=paste0(dir, output_file_prefix,"_", date, ".RDa"))

  #Write Files
  write.table(ASE.type, file=paste0(dir, output_file_prefix, "_", date, "_", 'event_types.txt'), sep = '\t', quote=FALSE, col.names=NA)
  write.table(FisherAnalyses, file=paste0(dir, output_file_prefix, "_", date, "_", 'FisherAnalyses.txt'), sep = '\t', quote=FALSE, col.names=NA)
  write.table(as.data.frame(junc.Outliers$TumorOverExpression), file = paste0(dir, output_file_prefix, "_", date, "_", 'TumorOverExpression.txt'), sep = '\t', quote=FALSE, col.names=NA)
  write.table(as.data.frame(junc.Outliers$TumorUnderExpression), file = paste0(dir, output_file_prefix, "_", date, "_",'TumorUnderExpression.txt'), sep = '\t', quote=FALSE, col.names=NA)
  write.table(splice_burden, file=paste0(dir, output_file_prefix, "_", date, "_", 'splice_burden.txt'), sep = '\t', quote=FALSE, col.names=NA)

  annotations_df <- annoGR2DF(geneAnnot)
  write.table(annotations_df, file = paste0(dir, output_file_prefix, "_", date, "_", 'gene_annotations.txt'), sep = '\t', quote = FALSE, col.names=NA)

  return(FisherAnalyses)
}
