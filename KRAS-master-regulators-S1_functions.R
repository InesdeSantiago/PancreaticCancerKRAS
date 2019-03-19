pkgTest <- function(x, repos="http://cran.rstudio.com/", tar=NULL)
{
  if (!is.null(repos)) {
      if (!require(x,character.only = TRUE)){install.packages(x,dep=TRUE,repos=repos)}
  }

  if (is.null(repos) & !is.null(tar)) {
    if (!require(x,character.only = TRUE)){
      install.packages(tar, repos=NULL, type="source")
      x <- tar
    }
  }

  if (is.null(repos) & is.null(tar)) {stop("cannot install tar file")}

  if(!require(x,character.only = TRUE))
    {
      source("http://bioconductor.org/biocLite.R")
      biocLite(x)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
}

########################################
#Convert Names: GeneSymbol_Entrez
########################################

#convertNames <- function(names2convert, mart) {
#    names <- mart$names[match(names2convert, mart$entrezgene)]
#    NAS <- which(is.na(names))
#    names[NAS] <- names2convert[NAS]
#    return(names)
#}


getentrez <- function(ES) {
  a <- sapply(ES,function (x) {r <- strsplit(x,"\\_")[[1]];
                                         if(length(r)==1) {return(r[[1]])}
                                         else{return(r[[ length(r) ]])}
                                           })
  names(a) <- NULL
  return(a)
}

getsymbol <- function(ES) {
  a <- sapply(ES,function (x) {r <- strsplit(x,"\\_")[[1]][[1]]})
  a <- a[a != ""]
  names(a) <- NULL
  return(a)
}


########################################
#Data preprocessing
########################################

plotQC_perstudy <- function(study_name,dataPLM,annot)
  {
  require(affy)
  cellfiles <- sampleNames(dataPLM)
  colvect <- annot[cellfiles,"tissue"]
  colvect[colvect=="PDAC"] <- "darkblue"
  colvect[colvect=="Normal"] <- "red"
  colvect[colvect=="cell_line"] <- "purple"

  Mbox(dataPLM,
    main = paste("RLE -", study_name),
    outline = FALSE, col = colvect, names=F,
    las = 3, whisklty = 0, staplelty = 0,ylim=c(-1,1))
  abline(h=0,lty=2)
  abline(h=0.15,lty=2,col="black")

  boxplot(dataPLM, main = paste("NUSE -", study_name),
          outline = FALSE, names=F, col = colvect,
          las = 3, whisklty = 0, staplelty = 0,ylim=c(0.9,1.3))
  abline(h=1,lty=2)
  abline(h=1.1,lty=2,col="black")

  #get outliers names
  RLE_stats <- RLE(dataPLM, type="stats")
  NUSE_stats <- NUSE(dataPLM, type="stats")
  RLE_bad <- colnames(RLE_stats)[which(as.numeric(RLE_stats["median",]) > 0.15)]
  Best_RLE <- colnames(RLE_stats)[order(abs(RLE_stats["median",]))[1:2]]
  NUSE_bad <- colnames(NUSE_stats)[which(as.numeric(NUSE_stats["median",]) > 1.1)]

  outliers <- list("RLE_bad" = RLE_bad, "NUSE_bad"=NUSE_bad, "Best_RLE"=Best_RLE, param=c("RLE"=0.15,"NUSE"=1.1))
  return(outliers)
}

#function to filter a gexp matrix
cv.filter <- function(gexp, ids){

  # compute CV
  meangx<-apply(gexp,1,mean,na.rm=TRUE)
  sdgx<-apply(gexp,1,sd,na.rm=TRUE)
  cvgx<-abs(sdgx/meangx)

  # aline ids
  ids <- data.frame(ids[match(rownames(gexp), ids[,1]),], "CV"=cvgx, stringsAsFactors=FALSE)

  # remove entries without annotation
  ids <- ids[!is.na(ids[,2]),]
  ids <- ids[ids[,2]!="",]
  print(paste("removed",nrow(gexp) - nrow(ids),"entries without annotation"))

  # remove duplicated genes by CV
  ids <- ids[sort.list(ids$CV,decreasing=TRUE),]
  unids<-unique(ids[,2])
  idx<-match(unids,ids[,2])
  print(paste("removed",nrow(ids) - length(idx),"duplicated entries"))
  ids<-ids[idx,]

  # update and return gexp matrix
  gexp <- gexp[rownames(gexp) %in% ids[,1],]
  rownames(gexp) <- ids[match(rownames(gexp), ids[,1]),2]
  ids <- ids[,-ncol(ids)]
  rownames(ids) <- ids[,2]
  ids <- ids[rownames(gexp),]
  return(list(gexp=gexp,ids=ids))
}

getGeneExpIDs <- function(gexp) {
  # Use the standard "db" file for the human U133plus2 array
  library(hgu133plus2hsentrezg.db)

  #Start with a list of all probesets on the array
  probesets.all = rownames(gexp)

  #Get some selected information
  gene.IDs = unlist(mget(probesets.all, hgu133plus2hsentrezgENTREZID, ifnotfound=NA))
  gene.symbols = unlist(mget(probesets.all, hgu133plus2hsentrezgSYMBOL, ifnotfound=NA))
  gene.names = unlist(mget(probesets.all, hgu133plus2hsentrezgGENENAME, ifnotfound=NA))

  #Print them all (and/or merge them with the expression matrix)
  probe.anno = data.frame(probesets.all, gene.IDs, gene.symbols, gene.names)
  return(probe.anno)
}


##############################################################################
#Functions to perform DEG - KRAS sig
##############################################################################

getfit <- function(copy, aWeights, targets, contrast, annotation) {
    #design matrix
    groups <- as.factor(targets$Sample_Group)
    design = model.matrix(~0 + groups)
    colnames(design) = sub("groups", "", colnames(design))

    #Estimate the correlation between measurements made on the same technical replicate:
    #Necessary to treat each technical replicated as a random effect
    corfit <- duplicateCorrelation(copy, design, block = colnames(copy))

    #This inter-technical replicate correlation is input into the linear model fit:
    fit = lmFit(copy, design, weights=aWeights, block = colnames(copy), cor=corfit$consensus.correlation)

    #Now we can make comparisons between the experimental conditions in the usual way:
    contrast.matrix = makeContrasts(contrasts=contrast,levels=design)

    #Then compute the contrast of interest and moderated t-test
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)

    #Add more columns to the table
    fit2$genes = annotation

    #sig
    sig <- topTable(fit2, adjust.method="BH", coef=1,number = Inf,resort.by="p")
    return(sig)

}

combineValues <- function(dframe, copy) {

  splitdata <- split(dframe[,c("IlluminaID","t")], as.character(dframe$hsapiens_entrezgene))
  splitdata <- lapply(splitdata, unique)

  #combine if more than 1 gene is found
  splitdata2 <- lapply( splitdata, function(d) {
    t <- d$t
    if (length(t) > 1) {t <- mean(t)}
    ps <- 2 *pt(abs(t), df=(ncol(copy)-1), lower.tail=F)
    res <- c("statistic"=t, "p.value"=ps)
    return(res)
  })
  splitdata2 <- do.call("rbind",splitdata2)
  rownames(splitdata2) <- names(splitdata)
  splitdata2
}

##############################################################################
#Functions to signature (VIPER format)
##############################################################################

getsighuman <- function(sig,annotTable,gexp) {

  copy <- gexp

  #merge
  sigM <- merge(unique(sig), unique(annotTable), by="IlluminaID", all.x=T)

  #Remove NA's and Badquality Probes
  nas <- which(is.na(sigM[,"hsapiens_entrezgene"]))
  sig2 <- sigM[-nas,] #DEG2 is signature without the NAS
  badquality <- which(sig2$PROBEQUALITY == "Bad")
  sig2 <- sig2[-badquality,] #delete bad quality

  #combine values when more than 1 human entrez gene exist
  v <- combineValues(sig2, copy)
  statistic <- matrix(v[,"statistic"])
  p.value <- matrix(v[,"p.value"])
  rownames(statistic) <- rownames(p.value) <- rownames(v)
  result <- list("statistic"=statistic, "p.value"=p.value)

  print (paste("initial dataset consisted of", length(unique(sigM$IlluminaID)), "unique probes"))
  print (paste("removed",  length(unique(sigM[nas,1])), "probes without any matching ID"))
  print (paste("removed", length(badquality), "lines of bad quality probes"))
  print (paste("final set consists of", length(unique(sig2$IlluminaID)), "unique probes"))

  print (paste("the signature contains", nrow(result[[1]]), "unique human ids"))

  return(result)
}

getZ <- function(sig,annotTable,copy) {

  #merge
  sigM <- merge(unique(sig), unique(annotTable), by="IlluminaID", all.x=T)

  #Remove NA's and Badquality Probes
  nas <- which(is.na(sigM[,"hsapiens_entrezgene"]))
  sig2 <- sigM[-nas,] #DEG2 is signature without the NAS
  badquality <- which(sig2$PROBEQUALITY == "Bad")
  sig2 <- sig2[-badquality,] #delete bad quality

  #combine values when more than 1 human entrez gene exist
  v <- combineValues(sig2, copy)
  statistic <- matrix(v[,"statistic"])
  rownames(statistic) <- rownames(v)

  return(statistic)
}

#some probes map to a NA and a Entrez HUman ID
#     IlluminaID          t   P.Value            ENSEMBL hsapiens_homolog_ensembl_gene hsapiens_entrezgene PROBEQUALITY
#16 ILMN_1212638 -0.9959256 0.3264195 ENSMUSG00000096319                          <NA>                  NA          Bad
#17 ILMN_1212638 -0.9959256 0.3264195 ENSMUSG00000062647               ENSG00000148303                6130          Bad

limmaNull <- function(copy, targets, annotTable, annotation, nPermutations=1000){

  #Estimates relative quality weights for each array in a multi-array experiment.
  groups <- as.factor(targets$Sample_Group)
  design = model.matrix(~0 + groups)
  colnames(design) = sub("groups", "", colnames(design))
  aWeights <- arrayWeights(copy,design=design)


  #Permutations
  Zmatrix <- lapply(1:nPermutations, function(i,annotTable) {
                   print(i)
                   set.seed(i)
                   #shuffle data
                   copy2 <- copy
                   perm <- sample(1:length(targets$Sample_Group))
                   aWeights2 <- aWeights[perm]
                   targets2 <- targets[perm,]

                   #get Z scores
                   sig <- getfit(copy2, aWeights2, targets2, contrast="Cre-Mock", annotation)
                   getZ(sig[,c("IlluminaID","t","P.Value")], annotTable,copy)


  }, annotTable=annotTable)
  Zmatrix <- do.call("cbind", Zmatrix)
  colnames(Zmatrix) <- 1:nPermutations
  return(Zmatrix)
}

##############################################################################
#Pathway analysis of KRASsig
##############################################################################

pathwayAnalysis <- function(fit, pvaluecutoff=0.01) {
    #prepare input data
    require(HTSanalyzeR)
    require(org.Mm.eg.db)
    require(GO.db)
    require(KEGG.db)
    require(snow)

    #First, create a named vector with the phenotype (t-statistic)
    #associated to each gene
    d <- na.omit(data.frame(fit[["ENTREZID"]], fit[["t"]]))
    data4enrich <- as.vector(d[,2])
    names(data4enrich) <- d[,1]

    #The hits as genes with adjusted p-value <= pvaluecutoff,
    geneNames <- fit[["ENTREZID"]]
    adjustedPvalues <- fit[["adj.P.Val"]]
    hits <- as.character(geneNames[adjustedPvalues <= pvaluecutoff])

    #define the gene set collections.
    GO_MF <- GOGeneSets(species="Mm", ontologies=c("MF"))
    GO_BP <- GOGeneSets(species="Mm", ontologies=c("BP"))
    PW_KEGG <- KeggGeneSets(species="Mm")
    ListGSC <- list(GO_MF=GO_MF, GO_BP=GO_BP, PW_KEGG=PW_KEGG)

    #add parallel computing to promote speed
    library(snow)
    options(cluster=makeCluster(4, "SOCK"))

    #overrepresentation and enrichment analyses
    gsca <- new("GSCA",
                  listOfGeneSetCollections=ListGSC,
                  geneList=data4enrich, hits=hits)
    gsca <- preprocess(gsca, species="Mm",
                  initialIDs="Entrez.gene",
                  keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                  orderAbsValue=FALSE, verbose=F)
    gsca <-analyze(gsca, para=list(pValueCutoff=pvaluecutoff,
                  pAdjustMethod="BH", nPermutations=100,
                  minGeneSetSize=15, exponent=1), verbose=F)


    if(is(getOption("cluster"), "cluster")) {
      stopCluster(getOption("cluster"))
      options(cluster=NULL)
    }
    return(gsca)
}

##############################################################################
#Preparing network data
##############################################################################

TFuse2regulon <- function(object, pvalueCutoff= 0.05, filterTFs=NULL, returnclassregulon=T) {

    if (class(object) == "TPC") {
        net <- object@results.pcor
        tfmode <- net$tfmode
        prob <- net$prob
        qval <- net$qval
    }
    if (class(object) == "list")    {
        net <- object
        if (!("tfmode" %in% names(net))) {stop("tfmode not found")}
        if (!("prob" %in% names(net))) {stop("prob not found")}
        if (!("qval" %in% names(net))) {stop("qval not found")}
        tfmode <- net$tfmode
        prob <- net$prob
        qval <- net$qval
    }

    if (is.null(filterTFs)) {TFs <- colnames(tfmode)}
    if (!is.null(filterTFs)) {TFs <- filterTFs; TFs <- intersect(TFs, colnames(tfmode))}
    scores <- prob[,TFs]
    mode1 <- tfmode[,TFs]

    aracne <- list()
    for (tf in TFs) {
        reg <- qval[,tf]
        which.reg <- reg <= pvalueCutoff
        likelihood <- scores[which.reg,match(tf, colnames(scores))]
        tfmode <- mode1[which.reg,match(tf, colnames(mode1))]
        aracne[[tf]] <- list("tfmode"=tfmode, "likelihood"=likelihood)
    }

    # removing missing data from the aracne regulon
    aracne <- aracne[names(aracne) != "NA"]
    aracne <- lapply(aracne, function(x) {
        filtro <- !(names(x$tfmode)=="NA" | is.na(x$tfmode) | is.na(x$likelihood))
        x$tfmode <- x$tfmode[filtro]
        x$likelihood <- x$likelihood[filtro]
        return(x)
    })
    regul <- aracne[sapply(aracne, function(x) length(names(x$tfmode)))>0]
    if (returnclassregulon) class(regul) <- "regulon"
    return(regul)
}


##############################################################################
#Compare Networks
##############################################################################


edgedistReturnFulldata <- function(TF, net1, net2, ...) {

    reg1 <- net1@results.pcor$z[,TF]
    reg2 <- net2@results.pcor$z[,TF]
    qval1 <- net1@results.pcor$qval[,TF]
    qval2 <- net2@results.pcor$qval[,TF]
    N <- intersect(names(reg1),names(reg2))
    reg1 <- reg1[N] #can only compare edges that exist in both gexp
    reg2 <- reg2[N]
    qval1 <- qval1[N]
    qval2 <- qval2[N]

    d <- cbind(reg1,reg2)
    rownames(d) <- names(reg1)
    return(d)
}

edgedistReturnSignEdges <- function(TF, net1, net2, combinedNet, pvaluecutoff=0.05, ...) {

    reg1 <- net1@results.pcor$z[,TF]
    reg2 <- net2@results.pcor$z[,TF]
    comb <- combinedNet$qval[,TF]
    sigedges <- names(comb)[comb <= pvaluecutoff]
    reg1 <- reg1[sigedges] #significant edges
    reg2 <- reg2[sigedges]
    d <- cbind(reg1,reg2)
    rownames(d) <- names(sigedges)
    return(d)
}


FETdist <- function(TF, net1, net2, pvaluecutoff=0.05) {
    reg1 <- net1@results.pcor$tn[,TF]
    qval1 <- net1@results.pcor$qval[,TF]
    bkg1 <- names(reg1)
    reg1 <- names(reg1)[qval1 <= pvaluecutoff]

    reg2 <- net2@results.pcor$tn[,TF]
    qval2 <- net2@results.pcor$qval[,TF]
    bkg2 <- names(reg2)
    reg2 <- names(reg2)[qval2 <= pvaluecutoff]

    bkg <- intersect(bkg1,bkg2)
    A <- factor(bkg %in% reg1, levels=c(TRUE, FALSE))
    B <- factor(bkg %in% reg2, levels=c(TRUE, FALSE))

    if (length(A)>15 & length(B) > 15) {
        tmp <- fisher.test(A, B, alternative="greater")$p.value
        return(-log10(tmp)) #return -log10 p-value
    }else{return(NA)}
}

getCor <- function(d, ...) {
    reg1 <- d[,1]
    reg2 <- d[,2]
    if (sd(reg1) != 0 & sd(reg2) != 0  & length(reg1) > 20 & length(reg2) > 20) {
        COR <- cor(reg1,reg2, ...)
        return(COR)
    }else{return(NA)}
}

compareNetworks <- function(net1, net2, combinedNetwork, pvaluecutoff = 0.05, ...) {

    TFtest <- intersect(colnames(net1@results.pcor$tn), colnames(net1@results.pcor$tn)) #TFs common in both networks
    TFtest <- intersect(TFtest, colnames(combinedNetwork$qval))

    #correlate edge Z scores per regulon - all edges
    d_all <- lapply(TFtest, edgedistReturnFulldata, net1 = net1, net2 = net2, method="pearson")
    names(d_all) <- TFtest
    CorPerReg_all <- lapply(d_all, function(x) {getCor(x, ...)})
    d_all <- do.call("rbind", d_all)
    CorOverall_all <- cor(d_all[,1], d_all[,2],...)

    #correlate edge Z scores per regulon - sig edges
    d_sig <- lapply(TFtest, edgedistReturnSignEdges, net1 = net1, net2 = net2, combinedNet=combinedNetwork, method="pearson", pvaluecutoff=pvaluecutoff)
    names(d_sig) <- TFtest
    CorPerReg_sig <- lapply(d_sig, function(x) {getCor(x, ...)})
    d_sig <- do.call("rbind", d_sig)
    CorOverall_sig <- cor(d_sig[,1], d_sig[,2],...)

    #fisher exact test to test regulon overlap - all edges
    FETPerReg <- unlist(lapply(TFtest, FETdist, net1 = net1, net2 = net2, pvaluecutoff=pvaluecutoff))
    names(FETPerReg) <- TFtest

    return(list("CorOverall_all"=CorOverall_all, "CorOverall_sig"=CorOverall_sig, "CorPerReg_all"=CorPerReg_all, "CorPerReg_sig"=CorPerReg_sig,"FETPerReg"=FETPerReg))

}

compareallNetworks <- function(allNetworks, combinedNetwork, pvaluecutoff = 0.05) {
    netnames <- names(allNetworks)
    compareRes <- list()
    for (i in 1:(length(netnames)-1)) {
      for (j in (i+1):length(netnames)) {
        if (i != j) {

          net1 <- allNetworks[[netnames[i]]]
          net2 <- allNetworks[[netnames[j]]]
          pairname <- paste0(netnames[i],".",netnames[j])
          print(pairname)
          a <- compareNetworks(net1, net2, combinedNetwork, pvaluecutoff = pvaluecutoff)
          compareRes[[pairname]] <- a
        }
      }
    }
    return(compareRes)
}


#plotting functions
do_x_range <- function(tmp, xlimit=NULL) {
  if (is.null(xlimit)) {x_range <- c(min(unlist(lapply(tmp, function(x) min(x$x)))),
                                  max(unlist(lapply(tmp, function(x) max(x$x)))))}else{x_range <- xlimit}
  return(x_range)
  }

#plot coefedges
plotscoresPerRegulon <- function(obj, type=c("density","boxplot"), ...) {

    colors <- rainbow(length(obj))

    if (type=="density") {
        #density plot showing distribution of the correlation coefficient

        tmp <- lapply(obj, function(x) density(na.omit(x),
                                                 from=range(na.omit(x))[1]))
        x_range <- do_x_range(tmp)
        yl <- max(unlist(lapply(tmp, function(x) max(x$y))))
        plot(tmp[[1]], col=colors[[1]], xlim=x_range,
            ylim=c(0, yl*1.05), lwd=2, ...)
        for (i in 2:length(tmp)) {lines(tmp[[i]], lwd=2, col=colors[[i]])}
        legend(x="topleft",legend=paste0(names(obj)), col=colors,
            lwd=4, cex=.8, lty=c(2,1))
    }

    if (type=="boxplot") {
        boxplot(obj, col=colors, las=2, ...)

    }

}


getActivity <- function(allNetworks, combinedNetwork, KRASsig, nullmodel, pvalueCutoff= 0.05) {
  require(viper)

  #individual networks
  nes <- list()
  nes.bt <- list()
  pv <- list()
  for (i in 1:length(allNetworks)) {
    netname <- names(allNetworks)[i]
    print(netname)
    net <- allNetworks[[i]]
    reguli <- TFuse2regulon(net, pvalueCutoff= pvalueCutoff, returnclassregulon=T)
    res <- msviper(KRASsig[["statistic"]], reguli, nullmodel, verbose=F)
    nes[[netname]] <- res$es$nes
    nes.bt[[netname]] <- res$es$nes.bt[,1]
    pv[[netname]] <- res$es$p.value
  }

  #combined
  reguli <- TFuse2regulon(combinedNetwork, pvalueCutoff= pvalueCutoff, returnclassregulon=T)
  res <- msviper(KRASsig[["statistic"]], reguli, nullmodel, verbose=F)
  combES <- res$es$nes
  combPV <- res$es$p.value

  #make a matrix - PVALUE RESULTS
  A <- lapply(pv, function(x) -log10(x))
  comb <- -log10(combPV)
  A <- c(A, "comb"=list(comb))
  allRegulators <- Reduce(union, lapply(A, names))
  a1 <- do.call("cbind", lapply(A, function(x) x[match(allRegulators,names(x))]))
  rownames(a1) <- allRegulators
  a1 <- a1[order(a1[,"comb"], decreasing=T),]

  #make a matrix - ES scroes
  A <- c(nes, "comb"=list(combES))
  allRegulators <- Reduce(union, lapply(A, names))
  a2 <- do.call("cbind", lapply(A, function(x) x[match(allRegulators,names(x))]))
  rownames(a2) <- allRegulators
  a2 <- a2[order(a2[,"comb"]),]

  return(list("pvalue.matrix"=a1, "NES.matrix"=a2))

}

##################################################
#Network visualiasation functions
##################################################
getRegulon <- function(netQ, TFs) {
    regulons <- lapply(TFs, function(tf) {
            reg <- netQ[,tf]
            return(names(reg)[reg <= 0.01])
    })
    names(regulons) <- TFs
    regulons
}

getsubsetofNetwork <- function(net, TFs=MRs) {
    require(reshape2)
    require(RedeR)

    #First, lets get the gene names of the targets of those MRs
    netQ <- net$qval
    regulons <- getRegulon(netQ, TFs=MRs)
    genes2draw <- unique(unlist(regulons, use.names=F))

    #now retrieve the network only for those genes/MRs edges
    netQ <- net$qval
    netQ <- netQ[genes2draw, names(regulons)]
    netQ <- melt(netQ)
    netQ <- netQ[netQ$value <= 0.01,]

    return(netQ)
}

plotcom <- function (net, MRs) {
    require(RedeR)
    require(igraph)

    subsetofNet <- getsubsetofNetwork(net, MRs)
    net4reder <- subsetofNet[,1:2]
    net4reder <- unique(net4reder)
    names(net4reder) <- c("target","TF")
    #net4reder$TF <- sapply(as.character(net4reder$TF), function(x) {strsplit(x, "_")[[1]][1]})
    #net4reder$target <- sapply(as.character(net4reder$target), function(x) {strsplit(x, "_")[[1]][1]})

    #comunity detection
    g <- graph.data.frame(net4reder, directed=F, vertices=NULL)
    g <- igraph:::simplify(g, remove.multiple = TRUE)
    fc <- fastgreedy.community(g)
    membership <- membership(fc)

    #three largest components
    a <- split(names(membership), membership)
    major <- as.numeric(names(sort(sapply(a, length),decreasing=T)[1:3]))

    #draw network
    library(RColorBrewer)
    library(reshape2)
    library(igraph)
    library(RedeR)

    rdp <- RedPort()
    calld(rdp, maxlag = 5000)
    resetd(rdp)
    g <- graph.data.frame(net4reder, directed=F, vertices=NULL)
    #addGraph(rdp, g, layout = NULL)
    vnames <- get.vertex.attribute(g, "name")
    V(g)$nodeFontSize <- ifelse(vnames %in% unique(net4reder$TF), 100,1)
    V(g)$nodeSize <- ifelse(vnames %in% unique(net4reder$TF), 100,30)
    V(g)$name <- sapply(as.character(V(g)$name), function(x) {ifelse (x %in% unique(net4reder$TF), strsplit(x, "_")[[1]][1], x)})

    colors <- membership[vnames]
    colors[!(colors %in% major)] <- "gray"
    colors[colors == major[1]] <- "deeppink"
    colors[colors == major[2]] <- "royalblue"
    colors[colors == major[3]] <- "springgreen"
    V(g)$nodeColor <- colors
    addGraph(rdp, g, layout = NULL)



    return(a[major])
}

#######################################################################
#Pathway analysis functions of the identified TF signatures
#######################################################################

Keggplot <- function(gene, universe, pvalueCutoff) {
  require(clusterProfiler)
  kk <- enrichKEGG(gene         = gene,
                 organism     = "human",
                 pAdjustMethod = "fdr",
                 pvalueCutoff = pvalueCutoff,
                 readable     = TRUE)
  return(kk)
}

reactomePAplot <- function(gene, universe, pvalueCutoff) {
  require(ReactomePA)
  kk <- enrichPathway(gene         = gene,
                 organism     = "human",
                 pAdjustMethod = "fdr",
                 pvalueCutoff = pvalueCutoff,
                 readable     = TRUE)
  return(kk)
}

barplotkk <- function(x, main) {
  x <- summary(x)
  x <- x[order(x$p.adjust, decreasing=T),]
  N=nrow(x)
  if (N > 15) {x <- x[ ((nrow(x)-14):nrow(x)),]}
  x$Description <- factor(x$Description, levels=x$Description)
  ggplot(x, aes(x=Description, y=-log10(p.adjust)))+
  geom_bar(stat="identity", color="black", fill="lightgray")+
  theme_minimal() + coord_flip() +
  xlab("") + ggtitle(main) + theme(axis.text.y = element_text(size=15),
  axis.text.x = element_text(size=14),  axis.title=element_text(size=14))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  require(ggplot2)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##############################################
#Run Viper
##############################################

runViper <- function(gexp, net, MR) {
  #regulMR
  regulMR <- TFuse2regulon(net, pvalueCutoff= 0.01, filterTFs=MR)
  require(viper)
  vpres <- viper(as.data.frame(gexp), regulMR)
  vpres
}

##########################################################
#test enrichment gene programes from Bailey et al., 2016
##########################################################

phypertest <- function(a1, a2, universe) {
  Overlap <- length(intersect(a1,a2))
  group2 <- length(a2)
  Totalnrgenes <- length(universe)
  group1 <- length(a1)
  phyper (Overlap, group2, Totalnrgenes-group2, group1, lower.tail=F)
}


enrch.GeneProgrames <- function(res, universe) {
  require(sigPathway)

  geneSymbols <- lapply(res, getsymbol)
  a <- gmxToG("../data/geneProgrames_Bailey.gmx", verbose = TRUE)
  GProgrames <- lapply(a,"[[", 3)
  names(GProgrames) <- sapply(a,"[[", 1)
  a <- data.frame("Hedgehog"= sapply(GProgrames,
                                     phypertest,
                                     a2 = geneSymbols[["Hedgehog"]],
                                     universe=universe),
           "Notch" = sapply(GProgrames,
                            phypertest,
                            a2 = geneSymbols[["Notch"]],
                            universe=universe),
           "CellCycle" = sapply(GProgrames,
                                phypertest,
                                a2 = geneSymbols[["CellCycle"]],
                                universe=universe))

  #correct pvalue
  a$GP <- rownames(a)
  a <- melt(a, id="GP")
  a$value <- -log10(p.adjust(a$value))
  a <- dcast(a, GP ~ variable)
  rownames(a) <- a$GP
  a <- a[,c(2,3,4)]

  #pheatmap
  annot <- data.frame("gene.programmes"=rownames(a))
  rownames(annot) <- rownames(a)
  ann_colors = rownames(a)
  names(ann_colors) <- ann_colors

  BREAKS <- c(0,seq(-log10(0.05), max(unlist(a)), length.out=10))
  COLOR <- colorRampPalette(rev(brewer.pal(n = 7,name ="Reds")))(length(BREAKS))
  pheatmap(a, scale="none",
           breaks=BREAKS,
           color = c("royalblue", rev(COLOR)),
           cluster_rows=T,
           annotation_row=annot,
           annotation_colors = list("gene.programmes"=ann_colors))
}


GeneProgramesPathways <- function(GProgrames, universe, fun) {
  require(org.Hs.eg.db)
  require(clusterProfiler)
  require(ReactomePA)
  x <- org.Hs.egSYMBOL2EG
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  entrezIDs <- lapply(GProgrames, function(gene) unlist(xx[gene], use.names=F))
  ck <- compareCluster(entrezIDs, fun=fun, organism="human")
  plot(ck, showCategory=10)
}

###########################
#Survival Analysis
###########################

mergewithGroups <- function(survmat, groups, groupLabels) {
	survmat <- survmat[(survmat[[1]] %in% names(groups)),]
	survmat <- cbind(survmat, "groups"=groups[match(survmat[[1]],names(groups))])

	survmat$groups <- as.character(survmat$groups)
	dic <- groupLabels
	names(dic) <- c("1","2","3")
	survmat$groups <- dic[survmat$groups]
	survmat$groups <- factor(survmat$groups,
                         levels=c("CellCycle", "Hedgehog", "Notch"))
	return(survmat)
}


########################################################################################
#Mut data - ICGC
########################################################################################

#create the table of OE per gene and p-value

getOE <- function(tableOfProp) {
    totalProp <- colSums(tableOfProp)
    proptions <- 100*tableOfProp["1",] / totalProp
    O <- proptions[c("CellCycle","Hedgehog","Notch")]
    E <- rowSums(tableOfProp[,c("CellCycle","Hedgehog","Notch")])
    E <- 100*E["1"]/sum(E)
    OE <- log2(O/E)
    return(OE)
}

getFreq <- function(tableOfProp) {
    totalProp <- colSums(tableOfProp)
    proptions <- 100*tableOfProp["1",] / totalProp
    O <- proptions[c("CellCycle","Hedgehog","Notch")]
    return(O)
}

getP <- function(tableOfProp) {
    totalProp <- rowSums(tableOfProp[,c("CellCycle","Hedgehog","Notch")])
    res <- c(fisher.test(cbind("group"=tableOfProp[,"CellCycle"], totalProp))$p.value,
           fisher.test(cbind("group"=tableOfProp[,"Hedgehog"], totalProp))$p.value,
           fisher.test(cbind("group"=tableOfProp[,"Notch"], totalProp))$p.value)
    #res <- c(fisher.test(tableOfProp[,c("CellCycle","Normal")])$p.value,
    #fisher.test(tableOfProp[,c("Hedgehog","Normal")])$p.value,
    #fisher.test(tableOfProp[,c("Notch","Normal")])$p.value)
    names(res) <- c("CellCycle","Hedgehog","Notch")
    return(res)
}

convertGeneSymbolsMutData <- function(longtable) {
  #convert geneSymbols
  genes <- longtable$gene
  data <- subset(mutdataGeneNames, hgnc_symbol != "" & !is.na(hgnc_symbol), select = c("ensembl_gene_id","hgnc_symbol"))
  geneSymbols <- data$hgnc_symbol[match(genes, data$ensembl_gene_id)]
  NAS <- which(is.na(geneSymbols))
  geneSymbols[NAS] <- genes[NAS]
  longtable$symbol <- geneSymbols
  longtable
}


applyTest <- function(mutdatabyTumour, tgroupsICGC) {

  tgroups <- factor(tgroupsICGC)
    levels(tgroups) <- c("CellCycle","Notch","Hedgehog")
    Tgroups <- as.character(tgroups)
    names(Tgroups) <- names(tgroups)
    #specimen <- read.delim("/Volumes/groups/Research/fmlab/public_folders/InesdeSantiago/Projects/PC_shivan/data_download_PDAC/Grimmond_clinical_data/specimen.PACA-AU.tsv", stringsAsFactors=F)
    #specimen <- specimen[,c("submitted_specimen_id", "specimen_type")]
    #normals1 <- specimen$submitted_specimen_id[specimen$specimen_type == "Normal - solid tissue"]
    #Ngroups <- rep("Normal", length(normals1))
    #names(Ngroups) <- normals1
    #groups <- c(Tgroups, Ngroups)
    groups <- Tgroups

  #test per gene
  #groups <- factor(groups, levels=c("CellCycle","Notch","Hedgehog","Normal"))
  allSamples <- names(groups)
  samplesNotPresent <- allSamples[!(allSamples %in% rownames(mutdatabyTumour))]
  samplesNotPresent <- matrix(0, nrow=length(samplesNotPresent), ncol=ncol(mutdatabyTumour))
  rownames(samplesNotPresent) <- allSamples[!(allSamples %in% rownames(mutdatabyTumour))]
  colnames(samplesNotPresent) <- colnames(mutdatabyTumour)
  gs <- colnames(mutdatabyTumour)
  mutdatabyTumour <- data.frame(rbind(mutdatabyTumour, samplesNotPresent), stringsAsFactors=F)
  colnames(mutdatabyTumour) <- gs
  groups <- groups[match(rownames(mutdatabyTumour), names(groups))]

  pb <- txtProgressBar(min = 0, max = ncol(mutdatabyTumour), style = 3)

  for (i in 1:ncol(mutdatabyTumour)) {
    geneName <- colnames(mutdatabyTumour)[i]
    x <- mutdatabyTumour[,i]
    x <- as.factor(x)
    tableOfProp <- table(x, groups)
    pvalues <- getP(tableOfProp)
    freq <- getFreq(tableOfProp)
    OE <- getOE(tableOfProp)
    res <- data.frame(freq, OE, pvalues)
    res$gene <- geneName
    res$group <- rownames(res)
    rownames(res) <- NULL
    if (i==1) {longtable <- res}else{longtable <- rbind(longtable, res)}

    #set progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)

  longtable$adj.pvalue <- p.adjust(longtable$pvalues)
  longtable$OE[is.na(longtable$OE)] <- min(longtable$OE[!(is.infinite(longtable$OE))], na.rm=T)
  longtable$OE[is.infinite(longtable$OE)] <- min(longtable$OE[!(is.infinite(longtable$OE))], na.rm=T)
  longtable <- convertGeneSymbolsMutData(longtable)
  longtable
}



