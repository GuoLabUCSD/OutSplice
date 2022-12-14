# Arguments:

# **data_file** - filepath to an R Data file containing output from the SplicingOutliers or TCGA_SplicingOutliers Functions.
#
# **NUMBER** - number of junctions to plot.  This can be top number of junctions (over or under expressed), or can be specific junctions in a list.  Default is 1.
#
# **junctions**- you can input the specific junction you want to graph (or vector of junctions), default is NULL
#
# **tail** - you can specify if you want top over or under expressed with tail='RIGHT' for junctions overexpressed in tumors, or tail='LEFT' for junctions underexpressed in tumors.  Default is NULL
#
# **p_value** - p-value threshold to use when plotting the top over or under-expressed junctions with "tail". Default = 0.05
#
# **GENE** - do you want to pick junctions based on a gene?  TRUE means you will pick all the junctions mapping to a certain gene.  FALSE means you do not pick based on the gene.
#
# **SYMBOL** - SYMBOL of gene you want to graph
#
# **makepdf** - do you want the graphs to be saved into a pdf?  Default is no (FALSE).
#
# **pdffile** - if you want to save to a pdf, you need to write the file path
#
# **tumcol** - color for tumors on graph. Default is red.
#
# **normcol** - color for normals on graph. Default is blue.

PlotJunctionData<-function(data_file, NUMBER=1, junctions=NULL, tail=NULL, p_value = 0.05, GENE=F, SYMBOL=NULL, makepdf=F, pdffile = NULL, tumcol='red', normcol='blue') {
  #library('NCBI2R')
  library('gplots')
  suppressPackageStartupMessages(load(data_file))
  ## if you want to make a pdf, this will be specified.  Stop/error if not specified
  if (makepdf==T) {
    if (is.null(pdffile)) stop('Need to specify PDF file path name. pdffile=...')
    pdf(pdffile)
  }
  ## Define which junctions are going to be graphed
  ## if you do not define which specific junctions you want to graph
  if (is.null(junctions)){
    if (GENE==F) {
      if (is.null(tail)){
        junctions<-row.names(junc.RPM)[1:NUMBER]
        #print(paste('Junction #',NUMBER, 'graphed.'))
      } else if (tail=='RIGHT'){
        topgenelist90.ogsa=FisherAnalyses[,'var2']
        toplist90<-FisherAnalyses[topgenelist90.ogsa,]
        if (sum(toplist90[,"FisherP2"]<p_value)>0){
          toplist90<-toplist90[toplist90[,"FisherP2"]<p_value,]
          toplist90<-toplist90[,"FisherP2"]
        }else{
          toplist90<-vector()
          print('no over expression outliers')
        }
        toplist90 <- as.data.frame(toplist90)
        junctions<-row.names(toplist90)[1:NUMBER]
        print(paste('Top', NUMBER, 'overexpressed junctions graphed.'))
      } else if (tail=='LEFT'){
        topgenelist10.ogsa=FisherAnalyses[,'var1']
        toplist10<-FisherAnalyses[topgenelist10.ogsa,]
        if (sum(toplist10[,"FisherP1"]<p_value)>0){
          toplist10<-toplist10[toplist10[,"FisherP1"]<p_value,]
          toplist10<-toplist10[,"FisherP1"]
        }else{
          toplist10<-vector()
          print('no under expression outliers')
        }
        toplist10 <- as.data.frame(toplist10)
        junctions<-row.names(toplist10)[1:NUMBER]
        print(paste('Top', NUMBER, 'underexpressed junctions graphed.'))
      }
    } else if (GENE==T){
      if (is.null(SYMBOL)) stop('Need to specify gene of interest by gene symbol. SYMBOL=\'...\'')
      junctions<-names(geneAnnot)[grep(SYMBOL, geneAnnot$SYMBOL)]
      print(paste(length(junctions), 'junction(s) for gene', SYMBOL, 'found.'))
    }
  }

  if (length(junctions)==0) stop('No junctions found.')
  count=0
  par(mfrow=c(2,2), oma=c(0,0,2,0))
  for (j in junctions)
  {
    ## first barplot is for raw junction data
    #print(junc.RPM.norm[j,])
    #View(junc.RPM.norm)
    print(j)
    phenotypes <- as.data.frame(pheno)
    samples <- order(phenotypes$pheno,junc.RPM.norm[j,],junc.RPM.norm[j,])
    #View(samples)
    barplot(log2(junc.RPM[j,samples]+1), cex.names=0.3, las=2, col = ifelse(phenotypes[samples,'pheno']=='Tumor',tumcol,normcol))
    title(sprintf('Junction Expression (log)\nskip? %s del? %s? ins? %s',
                  geneAnnot[j,]$skipping, geneAnnot[j,]$deletions,
                  geneAnnot[j,]$insertions),cex.main=0.7)

    ## this plot is for normalized junction expression - with RSEM.
    barplot(log2(junc.RPM.norm[j,samples]+1)-log2(NORM.RSEM.norm[j]+1),
            las=2, cex.names=0.3,
            col=ifelse(phenotypes[samples,'pheno']=='Tumor',tumcol,normcol))
    title('Junction expression \nNormalized by RSEM (log)', cex.main=0.7)

    ## this barplot is for RSEM of gene expression
    barplot(log2(RSEM[j,samples]+1),las=2, cex.names=0.3,
            col=ifelse(phenotypes[samples,'pheno']=='Tumor',tumcol,normcol))
    title(sprintf('%s \nRSEM gene expression',geneAnnot[j,]$SYMBOL), cex.main=0.7)

    plot.new()

    count=count+1
    #print(count)
  }



  if (makepdf==T) dev.off()
}
