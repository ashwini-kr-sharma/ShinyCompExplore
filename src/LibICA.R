###############################################################################
## Warp-ups for stable ICA: multi-run ICA providing consensus S and M matrices
## MIT License https://en.wikipedia.org/wiki/MIT_License
## (c) P.Nazarov, Luxembourg Institute of Health, 2018-06-30
## Thanks to T.Kaoma, M.Chepeleva and entire BIOMOD team of GenePro, LIH

## Supported by FNR CORE grant (C17/BM/11664971/DEMICS)
## ver.1.2009
###############################################################################

## ToDo:
## parallem functional annotation
## install reqired genomes 

if (!exists("runCrossValidation")) source("http://sablab.net/scripts/LibML.r")

## install packages required
{
print("----------------------------------------------------------------")
print("consICA: Checking required packages and installing missing ones")
packages = c("fastICA","sm","vioplot","pheatmap","gplots","foreach","survival","topGO","org.Hs.eg.db")
if(.Platform$OS.type=="unix")  packages = c(packages,"doMC")
if(.Platform$OS.type=="windows")  packages = c(packages,"doSNOW")
print("Needed:")
print(packages)
print("Absent:")
print(packages[!packages %in% rownames(installed.packages())])

if (!requireNamespace("BiocManager", quietly = TRUE)) 
	install.packages("BiocManager")

for (pk in packages){
	if (! pk %in% rownames(installed.packages()))
		BiocManager::install(pk)
}
rm(packages,pk)
}
##=============================================================================
## runICA - runs consensus ICA with parallelization. 
## Estimates S and M matrices such that X=SxM, S has statistically independent columns
##	X      - data matrix (features - rows, samples - columns)
##	ncomp  - number of components
##	ntry   - number of runs
##	show.every - how often do we show the progress messages (disabled for ncores>1)
##	ncores - number of cores to be set for parallel calculation
##	exclude - are samples excluded during multiple run?
##	Returns: list with 
##			X - data
##			S,M - consensus metagene and weight matrix
##			mse - error ||X-SxM||
##			mr2 - mean R2 between rows of M
##			stab - stability, mean R2 between consistent columns of S in multiple tries
##-----------------------------------------------------------------------------
runICA = function(
	X,
	ncomp=3, 
	ntry = 1, 
	show.every=1, 
	filter.thr = NULL,
	ncores=1, 
	reduced=FALSE, 
	savecomp=NULL, 
	exclude=TRUE, 
	fun="logcosh", 
	alg.typ="deflation",
	method = c("R","C"),
	seed=NULL)
{
	if (!is.null(seed)) set.seed(seed)
	
	if (ncores > ntry) ncores = ntry

		
	## install packages if absent
	if (!"fastICA" %in% rownames(installed.packages())){
		print("Cannot find `fastICA` - installing it from Bioconductor")
		source("https://bioconductor.org/biocLite.R")
		biocLite("fastICA")
	}
	if (!"sm" %in% rownames(installed.packages())){
	  print("Cannot find `fastICA` - installing it from Bioconductor")
	  source("https://bioconductor.org/biocLite.R")
	  biocLite("sm")
	}
	if (ncores>1 &.Platform$OS.type=="unix" & !"doMC" %in% rownames(installed.packages())){
		print("Cannot find `doMC` - installing it from Bioconductor")
		source("https://bioconductor.org/biocLite.R")
		biocLite("doMC")
	}
	if (ncores>1 &.Platform$OS.type=="windows" & !"doSNOW" %in% rownames(installed.packages())){
		print("Cannot find `doSNOW` - installing it from Bioconductor")
		source("https://bioconductor.org/biocLite.R")
		biocLite("doSNOW")
	}
		
	require(fastICA)
	## create containers
	X = as.matrix(X)
	if (!is.null(filter.thr)) X = X[apply(X,1,max)>filter.thr,]
	Res = list()
	Res$X = X
	S = list()
	M = list()
	## metagenes
	S[[1]] = matrix(nrow=nrow(X),ncol=ncomp)
	rownames(S[[1]]) = rownames(X)
	colnames(S[[1]]) = sprintf("ic.%d",1:ncomp)
	## mixing matrix
	M[[1]] = matrix(nrow=ncomp,ncol = ncol(X))
	colnames(M[[1]]) = colnames(X)
	rownames(M[[1]]) = sprintf("ic.%d",1:ncomp)
	itry=1
	Res$S = S[[1]]
	Res$M = M[[1]]
	Res$mse = NA  ## mean square error bw X and S*M
	Res$mr2 = NA  ## mean correlation bw mixing profiles
	Res$n3 = NA  ## mean number of elements in |S| over 3
	## do multiple tries
	idx.excludeSamples = sample(1:ncol(X),ntry, replace = (ntry > ncol(X)))
	if (ntry==1 | (!exclude) ) idx.excludeSamples = integer(0)
	itry = 1
	
	###############################
	## Parallel section starts
	##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	if (ncores > 1) {
		require(foreach)
		if(.Platform$OS.type=="unix") {
			require(doMC)
			registerDoMC(ncores)
		}
		if(.Platform$OS.type=="windows") {
			require(doSNOW)
			cl = makeCluster(ncores) 
			registerDoSNOW(cl)
		}
	}
	
	cat("*** Starting",ifelse(ncores>1,"parallel",""),"calculation on",ncores,"core(s)...\n")
	cat("*** System:",.Platform$OS.type,"\n")
	cat("***",ncomp,"components,",ntry,"runs,",nrow(X),"features,",ncol(X),"samples.\n")
	cat("*** Start time:",as.character(Sys.time()),"\n")
	t0=Sys.time()
	flush.console()
	## multi run ICA
	if (ncores > 1) {
		MRICA = foreach(itry=1:ntry) %dopar% {
			require(fastICA)
			SP = Res$S + NA
			MP = Res$M + NA
			if (length(idx.excludeSamples) == 0){
				ic = fastICA(X, n.comp = ncomp, alg.typ=alg.typ ,fun=fun,method=method)
				SP[,] = ic$S
				MP[,] = ic$A
			}else{
				x = X[,-idx.excludeSamples[itry]]
				ic = fastICA(x, n.comp = ncomp, alg.typ=alg.typ ,fun=fun,method=method)
				SP[,]= ic$S
				MP[,-idx.excludeSamples[itry]] = ic$A
				MP[is.na(MP)] = 0
			}
			gc()
			return(list(S=SP,M=MP))
		}
	}else{
		cat("Execute one-core analysis, showing progress every",show.every,"run(s)\n")
		require(fastICA)
		MRICA = list()
		for(itry in 1:ntry){
			MRICA[[itry]] = list()
			MRICA[[itry]]$S = Res$S + NA
			MRICA[[itry]]$M = Res$M + NA
			if (length(idx.excludeSamples) == 0){
				ic = fastICA(X, n.comp = ncomp, alg.typ=alg.typ ,fun=fun,method=method)
				MRICA[[itry]]$S[,] = ic$S
				MRICA[[itry]]$M[,] = ic$A
			}else{
				x = X[,-idx.excludeSamples[itry]]
				ic = fastICA(x, n.comp = ncomp, alg.typ=alg.typ ,fun=fun)
				MRICA[[itry]]$S[,]= ic$S
				MRICA[[itry]]$M[,-idx.excludeSamples[itry]] = ic$A
				MRICA[[itry]]$M[is.na(MRICA[[itry]]$M)] = 0
			}
			if (itry%%show.every == 0) {
				cat("try #",itry,"of",ntry,"\n")
				flush.console()
			}
		}
	}
	
	if(.Platform$OS.type=="windows" & ncores>1)  stopCluster(cl)
	
	cat("*** Done!", "\n")
	cat("*** End time:",as.character(Sys.time()),"\n")
	print(Sys.time()-t0)
	
	##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	## ...end of parallel section
	###############################
	
	S = lapply(MRICA,function(x)x$S)
	M = lapply(MRICA,function(x)x$M)
	rm(MRICA)
	
	#print(str(S))
	
	## this is for debugging reasons - save partial results
	if (!is.null(savecomp)) save(list=ls(),file=savecomp)
	
	cat("Calculate ||X-SxM|| and r2 between component weights\n")
	for (itry in 1:ntry){
		Res$mse[itry] = sum((X - mean(X) - S[[itry]]%*%M[[itry]])^2)/ ncol(X) / nrow(X)
		## modified 2019-10-07/ Check!!
		Res$mr2[itry] = mean((cor(t(M[[itry]]))^2)[upper.tri(matrix(0,nrow=ncomp,ncol=ncomp))])
		Res$n3[itry] = sum(apply(abs(S[[itry]])>3,2,sum))/ncomp
	}
	## who is the best: min correlated (hoping that we capture majority of independent signals
	Res$i.best = which.min(Res$mr2)
	## simply smallest error? #Res$i.best = which.min(Res$mse) - no!
	if (length(Res$i.best) == 0 ) Res$i.best = 1
	## if only one try - return
	if (ntry == 1) {
		Res$S = S[[1]]
		Res$M = M[[1]]
		return(Res)
	}
	
	## correlate results
	## s.cor - to which ic of the BEST decomposition we should address? 
	cat("Correlate rows of S between tries\n")
	s.cor = matrix(nrow=ntry,ncol=ncomp)
	s.cor[Res$i.best,] = 1:ncomp
	itry = 1
	for (itry in (1:ntry)[-Res$i.best]) {
		r = cor(S[[itry]],S[[Res$i.best]])
		s.cor[itry,] = apply((r)^2,2,which.max)
		for (ic in 1:ncomp)
			s.cor[itry,ic] = s.cor[itry,ic] * sign(r[s.cor[itry,ic],ic])
	}
	## build consensus S, M
	
	cat("Build consensus ICA\n")
	Res$S[,] = S[[1]]
	Res$M[,] = M[[1]]
	itry=2
	for (itry in 2:ntry){ ## itry=2, because 1 is already there
		for (ic in 1:ncomp) {
			Res$S[,ic] = Res$S[,ic] + S[[itry]][,abs(s.cor[itry,ic])]* sign(s.cor[itry,ic])
			Res$M[ic,] = Res$M[ic,] + M[[itry]][abs(s.cor[itry,ic]),]* sign(s.cor[itry,ic])
		}
	}
	Res$S = Res$S / ntry
	Res$M = Res$M / ntry
	
	## use consensus S, M to analyze stability
	cat("Analyse stability\n")
	Res$stab = s.cor + NA
	for (itry in (1:ntry)) {
		Res$stab[itry,] = diag(cor(Res$S,S[[itry]][,abs(s.cor[itry,])])^2)
	}
	colnames(Res$stab) = colnames(Res$S) 
	
	
	if (reduced) {
		Res$X = NULL
		Res$mse = NULL
		Res$n3 = NULL
		Res$i.best = NULL
	}
	
	gc()
	
	return(Res)
}


runICAonce = function(X,ncomp=3, filter.thr = NULL,reduced=FALSE, fun="logcosh", alg.typ="deflation"){ 
	## install packages if absent
	if (!"fastICA" %in% rownames(installed.packages())){
		print("Cannot find `fastICA` - installing it from Bioconductor")
		source("https://bioconductor.org/biocLite.R")
		biocLite("fastICA")
	}
	require(fastICA)
	## create containers
	X = as.matrix(X)
	if (!is.null(filter.thr)) X = X[apply(X,1,max)>filter.thr,]
	Res = list()
	if (!reduced) Res$X = X
	Res$S = matrix(nrow=nrow(X),ncol=ncomp)
	rownames(Res$S) = rownames(X)
	colnames(Res$S) = sprintf("ic.%d",1:ncomp)
	## mixing matrix
	Res$M = matrix(nrow=ncomp,ncol = ncol(X))
	colnames(Res$M) = colnames(Res$X)
	rownames(Res$M) = sprintf("ic.%d",1:ncomp)
	ic = fastICA(X, n.comp = ncomp, alg.typ =alg.typ, fun=fun)
	Res$S[,] = ic$S
	Res$M[,] = ic$A
	colanmes(Res$M) = colnames(X)
	return(Res)
}

setOrientation = function(IC, GO){
	Res = list()
	Res$IC = IC
	Res$GO = GO
	Res$compsign = double(ncol(IC$S))
	for(ic in 1:ncol(IC$S)){
		pos = min(GO$GOBP[[ic]]$pos$FDR)
		neg = min(GO$GOBP[[ic]]$neg$FDR)
		if(neg < pos) {
			cat("rotating IC",ic,"\n")
			Res$compsign[ic] = -1
			Res$IC$S[,ic] = -Res$IC$S[,ic]
			Res$IC$M[ic,] = -Res$IC$M[ic,]
			for (db in c("GOBP","GOCC","GOMF")){
				tmp = Res$GO[[db]][[ic]]$pos
				Res$GO[[db]][[ic]]$pos = Res$GO[[db]][[ic]]$neg
				Res$GO[[db]][[ic]]$neg = tmp
			}
		} else {
			Res$compsign[ic] = 1
		}
	}
	return(Res)
}

##=============================================================================
## testing ICA: for debug and investigation only
##=============================================================================
# testParallelICA=function(data.file = "H:/BioinformaticsGroup/Members/Peter/astro.RData") {
	# load(data.file)#load("B:/Data/astro.RData")
	# ncomp=8
	# ntry = 1000
	# show.every=10
	# filter.thr = 7
	# ncores = 1
	# IC = runICA(X,ncomp=ncomp,ntry=ntry,show.every=show.every,filter.thr = 7,ncores=1) #6.45
	# IC = runICA(X,ncomp=ncomp,ntry=ntry,show.every=show.every,filter.thr = 7,ncores=2) #3.58
	# IC = runICA(X,ncomp=ncomp,ntry=ntry,show.every=show.every,filter.thr = 7,ncores=3) #3.03
	# IC = runICA(X,ncomp=ncomp,ntry=ntry,show.every=show.every,filter.thr = 7,ncores=4) #2.75
	# IC = runICA(X,ncomp=ncomp,ntry=ntry,show.every=show.every,filter.thr = 7,ncores=5) #2.73
	# IC = runICA(X,ncomp=ncomp,ntry=ntry,show.every=show.every,filter.thr = 7,ncores=6) #2.78
	# IC = runICA(X,ncomp=ncomp,ntry=ntry,show.every=show.every,filter.thr = 7,ncores=8) #2.99
	# IC = runICA(X,ncomp=ncomp,ntry=ntry,show.every=show.every,filter.thr = 7,ncores=10) #2.98
# }

getTopIdx=function(x,n){
	return(order(x,na.last=TRUE,decreasing=TRUE)[1:n])
}

##=============================================================================
## getGenesICA - Get influential genes after ICA
##	IC - IC list from runICA
##	alpha - FDR level for top-contibuting gene selection
##-----------------------------------------------------------------------------
getGenesICA = function(IC,alpha = 0.05){
	Genes = list()
	icomp=1
	for (icomp in 1:ncol(IC$S)){
		z = IC$S[,icomp]
		z = (z - median(z)) / mad(z)
		pv = pnorm(z)
		pv[pv>0.5] = 1-pv[pv>0.5]
		fdr=p.adjust(pv,method="BH")
		fdr.pos = fdr
		fdr.neg = fdr
		fdr.pos[fdr>alpha | z<0]=1
		fdr.neg[fdr>alpha | z>0]=1
		genes = rownames(IC$S) ## ensembl
		Genes[[sprintf("ic%02d",icomp)]]=list()
		Genes[[sprintf("ic%02d",icomp)]]$pos = data.frame(genes = genes[fdr.pos<1],fdr = fdr[fdr.pos<1], stringsAsFactors=FALSE)
		Genes[[sprintf("ic%02d",icomp)]]$neg = data.frame(genes = genes[fdr.neg<1],fdr = fdr[fdr.neg<1], stringsAsFactors=FALSE)
	}
	return(Genes)
}

##=============================================================================
## getGO - assigns IC signatures to Gene Ontologies
##	IC - IC structure from runICA
##	alpha - FDR level (both for top-contributing genes and GO enrichment)
##	lib.url - place where enrichGO.r is located (topGO warp-up)
##-----------------------------------------------------------------------------
getGO = function(IC,alpha = 0.05, genenames=NULL, lib.url = "http://sablab.net/scripts",genome="org.Hs.eg.db"){
	## get enrichGO library
	if (!exists("enrichGO")) source(paste(lib.url,"enrichGO.r",sep="/"))
	## structures
	Genes = list()
	GOBP = list()
	GOCC = list()
	GOMF = list()
	#REPA = list()
	## main loop
	icomp=1
	for (icomp in 1:ncol(IC$S)){
		pv = pnorm(IC$S[,icomp],m=median(IC$S[,icomp]),s=mad(IC$S[,icomp]))
		pv[pv>0.5] = 1-pv[pv>0.5]
		fdr=p.adjust(pv,method="BH")
		fdr.pos = fdr
		fdr.neg = fdr
		fdr.pos[fdr>alpha | scale(IC$S[,icomp])<0]=1
		fdr.neg[fdr>alpha | scale(IC$S[,icomp])>0]=1
		
		if (is.null(genenames)) {
			genes = rownames(IC$S) 
		}else{
			if(!is.null(names(genenames))){
				genes = genenames[rownames(IC$S)]
			}else{
				genes = genenames
			}
		}
		
		GOBP[[sprintf("ic%02d",icomp)]] = list()
		GOCC[[sprintf("ic%02d",icomp)]] = list()
		GOMF[[sprintf("ic%02d",icomp)]] = list()
		#REPA[[sprintf("ic%02d",icomp)]] = list()
		cat("=========================================================\nComponent",icomp,"\n\n")
		# tst = enrichPathway(gene = genes[fdr.pos<1],universe =genes, qvalueCutoff=0.01, readable=T)
		# if (!is.null(tst)){ 
			# REPA[[icomp]]$pos = tst@result[,c("Description","qvalue")]
			# names(REPA[[icomp]]$pos) = c("Term","FDR")
		# }else{REPA[[icomp]]$pos = ""}
		# tst = enrichPathway(gene = en.genes[fdr.neg<1],universe = en.genes, qvalueCutoff=0.01, readable=T)
		# if (!is.null(tst)){ 
			# REPA[[icomp]]$neg = tst@result[,c("Description","qvalue")]
			# names(REPA[[icomp]]$pos) = c("Term","FDR")
		# }else{REPA[[icomp]]$neg = ""}
		GOBP[[icomp]]$pos = enrichGO(genes = genes,fdr = fdr.pos,thr.fdr=alpha,db="BP",id= c("entrez", "ensembl", "symbol", "genename", "unigene"),genome=genome)
		GOCC[[icomp]]$pos = enrichGO(genes = genes,fdr = fdr.pos,thr.fdr=alpha,db="CC",id= c("entrez", "ensembl", "symbol", "genename", "unigene"),genome=genome)
		GOMF[[icomp]]$pos = enrichGO(genes = genes,fdr = fdr.pos,thr.fdr=alpha,db="MF",id= c("entrez", "ensembl", "symbol", "genename", "unigene"),genome=genome)
		cat("\n\nBP:",sum(GOBP[[icomp]]$pos$FDR<alpha),"enriched\n");#print(GOBP[[icomp]]$pos[GOBP[[icomp]]$pos$FDR<0.01,])
		cat("\n\nCC:",sum(GOCC[[icomp]]$pos$FDR<alpha),"enriched\n");#print(GOCC[[icomp]]$pos[GOCC[[icomp]]$pos$FDR<0.01,])
		cat("\n\nMF:",sum(GOMF[[icomp]]$pos$FDR<alpha),"enriched\n");#print(GOMF[[icomp]]$pos[GOMF[[icomp]]$pos$FDR<0.01,])
		GOBP[[icomp]]$neg = enrichGO(genes = genes,fdr = fdr.neg,thr.fdr=alpha,db="BP",id= c("entrez", "ensembl", "symbol", "genename", "unigene"),genome=genome)
		GOCC[[icomp]]$neg = enrichGO(genes = genes,fdr = fdr.neg,thr.fdr=alpha,db="CC",id= c("entrez", "ensembl", "symbol", "genename", "unigene"),genome=genome)
		GOMF[[icomp]]$neg = enrichGO(genes = genes,fdr = fdr.neg,thr.fdr=alpha,db="MF",id= c("entrez", "ensembl", "symbol", "genename", "unigene"),genome=genome)
		cat("\n\nBP:",sum(GOBP[[icomp]]$neg$FDR<alpha),"enriched\n");#print(GOBP[[icomp]]$neg[GOBP[[icomp]]$neg$FDR<0.01,])
		cat("\n\nCC:",sum(GOCC[[icomp]]$neg$FDR<alpha),"enriched\n");#print(GOCC[[icomp]]$neg[GOCC[[icomp]]$neg$FDR<0.01,])
		cat("\n\nMF:",sum(GOMF[[icomp]]$neg$FDR<alpha),"enriched\n");#print(GOMF[[icomp]]$neg[GOMF[[icomp]]$neg$FDR<0.01,])
		flush.console()
	}
	GO = list(GOBP = GOBP, GOCC = GOCC, GOMF=GOMF)#, REPA = REPA)
	#if (reCalcICA)	save(list=c("X","Meta","Var","ncomp","IC","Genes","FEA"),file=sprintf("ICA_nc=%02d_%s_GO.RData",ncomp,measure))
	return(GO)
}

##=============================================================================
## saveGO - saves GO and Genes into text files
##-----------------------------------------------------------------------------
saveGO = function(Genes,GO,folder=""){
	if(!exists("sortDataFrame")) source("http://sablab.net/scripts/sortDataFrame.r")
	ncomp = length(Genes)
	if (folder=="") folder = sprintf("ica_%02d",ncomp)
	dir.create(folder,showWarnings=FALSE)
	icomp=1
	for (icomp in 1:ncomp){
		if (nrow(Genes[[icomp]]$pos)>0)
			tab = sortDataFrame(Genes[[icomp]]$pos,"fdr")
		write.table(tab,file=sprintf("%s/ic%03d_genes_pos.txt",folder,icomp), row.names=FALSE, quote = FALSE, sep="\t")
		if (nrow(Genes[[icomp]]$neg)>0)
			tab = sortDataFrame(Genes[[icomp]]$neg,"fdr")
		write.table(tab,file=sprintf("%s/ic%03d_genes_neg.txt",folder,icomp), row.names=FALSE, quote = FALSE, sep="\t")
		goname="GOBP";d="pos";
		for (goname in c(names(GO))){
			for (d in c("pos","neg")){
				igo = GO[[goname]][[icomp]][[d]]$FDR<0.01
				fl = sprintf("%s/ic%03d_%s_%s.txt",folder,icomp,goname,d)
				write.table(GO[[goname]][[icomp]][[d]][igo,c("GO.ID","Term","FDR")],file = fl,sep="\t",row.names=FALSE,quote=FALSE)
			}
		}
	}
}
##=============================================================================
## saveGenes - saves Genes into text files
##-----------------------------------------------------------------------------
saveGenes = function(Genes,folder=""){
	if(!exists("sortDataFrame")) source("http://sablab.net/scripts/sortDataFrame.r")
	ncomp = length(Genes)
	if (folder=="") folder = sprintf("ica_%02d",ncomp)
	dir.create(folder,showWarnings=FALSE)
	icomp=1
	for (icomp in 1:ncomp){
		tab = sortDataFrame(Genes[[icomp]]$pos,"fdr")
		write.table(tab,file=sprintf("%s/ic%03d_genes_pos.txt",folder,icomp), row.names=FALSE, quote = FALSE, sep="\t")
		tab = sortDataFrame(Genes[[icomp]]$neg,"fdr")
		write.table(tab,file=sprintf("%s/ic%03d_genes_neg.txt",folder,icomp), row.names=FALSE, quote = FALSE, sep="\t")
	}
}

##=============================================================================
## saveReport - saves PDF report on each IC
##-----------------------------------------------------------------------------
saveReport = function(IC, Genes=NULL, GO=NULL, Var=NULL, surv=NULL, genenames=NULL,file = sprintf("report_ICA_%d.pdf",ncol(IC$S)), 
                      main = "Component # %d (stability = %.3f)", show.components = 1:ncol(IC$S)){
	require(pheatmap)
	require(gplots)
	if(!exists("drawTable"))  source("http://sablab.net/scripts/drawTable.r")
	if(!exists("violinplot")) source("http://sablab.net/scripts/violinplot.r")
	if(is.null(Genes)) Genes = getGenesICA(IC,alpha=0.01)
	if(!is.null(Var)) if (class(Var) != "data.frame") Var = as.data.frame(Var)
	
	try(dev.off(),silent=TRUE)
	try(dev.off(),silent=TRUE)
	try(dev.off(),silent=TRUE)
	
	ncomp = ncol(IC$S)
	pdf(file,width=8.3, height=11.7,onefile=TRUE)
	icomp=1
	for (icomp in 1:ncomp){
		for (direct in c("neg","pos"))
			if (nrow(Genes[[icomp]][[direct]])>1)
				Genes[[icomp]][[direct]] = sortDataFrame(Genes[[icomp]][[direct]],"fdr")
	}
	for (icomp in show.components){
		cat("Working with component #",icomp,"\n")
		par(mfcol=c(1,1),mar=c(3,3,2,1))
		plot.new()
		title(sprintf(main,icomp,mean(IC$stab[,icomp])),cex.main=0.8)
		#title(bquote("Component #" ~ .(icomp) ~ "( mean"~ R^2 ~ "="~ .(sprintf("%.3f",mean(IC$stab[,icomp])))~")" ),cex.main=0.8)
		## gene signature
		par(fig=c(0,0.2,0.85,1),new=TRUE,mar=c(2,2,2,1))
		plot(sort(IC$S[,icomp]),col = "#0000FF", type="l", ylab=("involvement"),xlab="genes",las=2,cex.axis=0.4, main="Metagene\n(involvement of features)",cex.main=0.6)
		text(0,0,"negative",adj=c(0,-1),cex=0.5,col="#000088")
		text(nrow(IC$S),0,"positive",adj=c(1,1),cex=0.5,col="#880000")
		par(fig=c(0,0.24,0,0.85),new=TRUE,mar=c(2,1,0,1))
		plot.new()
		text(0,1,paste(nrow(Genes[[icomp]]$neg),"\nnegative"),font=2,adj=c(0,0),cex=0.5,col="#000088")
		if (nrow(Genes[[icomp]]$neg)>0){
			txt = Genes[[icomp]]$neg$genes
			if (!is.null(genenames)) txt = genenames[txt]
			for (i in 1:nrow(Genes[[icomp]]$neg)) text(0,0.99-i/80,txt[i],col="#000088",adj=c(0,0),cex=0.4)
		}
		rect(0.5,0,1,1,col="white",border=NA)
		text(0.5,1,paste(nrow(Genes[[icomp]]$pos),"\npositive"),font=2,adj=c(0,0),cex=0.5,col="#880000")
		if (nrow(Genes[[icomp]]$pos)>0){
			#for (i in 1:nrow(Genes[[icomp]]$pos)) text(0.5,1-i/60,Genes[[icomp]]$pos$genes[i],col="#880000",adj=c(0,0),cex=0.5)
			txt = Genes[[icomp]]$pos$genes
			if (!is.null(genenames)) txt = genenames[txt]
			for (i in 1:nrow(Genes[[icomp]]$pos)) text(0.5,0.99-i/80,txt[i],col="#880000",adj=c(0,0),cex=0.4)
		}
		
		if (!is.null(GO)){
			par(fig=c(0.2,0.6,0.7,1),new=TRUE,mar=c(2,2,2,0));plot.new()
			tab = GO$GOBP[[icomp]]$neg[,c("Term","FDR")]
			tab = tab[tab$FDR<0.1,]
			tab$FDR =sprintf("%.2e",tab$FDR)
			text(0,1,sprintf("GO:BP neg : %d terms(FDR<0.01)",nrow(tab)),col="#000088",font=2,adj=c(0,0),cex=0.6)
			if (nrow(tab)>0) drawTable(tab[1:min(20,nrow(tab)),],x0=0,y0=0.98,dx=c(0.8,0.2),dy=0.04,row.names=FALSE,cex=0.5,col="#000088")
			par(fig=c(0.6,1,0.7,1),new=TRUE,mar=c(2,2,2,0));plot.new()
			tab = GO$GOBP[[icomp]]$pos[,c("Term","FDR")]
			tab = tab[tab$FDR<0.1,]
			tab$FDR =sprintf("%.2e",tab$FDR)
			text(0,1,sprintf("GO:BP pos : %d terms(FDR<0.01)",nrow(tab)),col="#880000",font=2,adj=c(0,0),cex=0.6)
			if (nrow(tab)>0) drawTable(tab[1:min(20,nrow(tab)),],x0=0,y0=0.98,dx=c(0.8,0.2),dy=0.04,row.names=FALSE,cex=0.5,col="#880000")
			
			par(fig=c(0.2,0.6,0.506,0.806),new=TRUE,mar=c(2,2,2,0));plot.new()
			tab = GO$GOCC[[icomp]]$neg[,c("Term","FDR")]
			tab = tab[tab$FDR<0.1,]
			tab$FDR =sprintf("%.2e",tab$FDR)
			text(0,1,sprintf("GO:CC neg : %d terms(FDR<0.01)",nrow(tab)),col="#000088",font=2,adj=c(0,0),cex=0.6)
			if (nrow(tab)>0) drawTable(tab[1:min(10,nrow(tab)),],x0=0,y0=0.98,dx=c(0.8,0.2),dy=0.04,row.names=FALSE,cex=0.5,col="#000088")
			par(fig=c(0.6,1,0.506,0.806),new=TRUE,mar=c(2,2,2,0));plot.new()
			tab = GO$GOCC[[icomp]]$pos[,c("Term","FDR")]
			tab = tab[tab$FDR<0.1,]
			tab$FDR =sprintf("%.2e",tab$FDR)
			text(0,1,sprintf("GO:CC pos : %d terms(FDR<0.01)",nrow(tab)),col="#880000",font=2,adj=c(0,0),cex=0.6)
			if (nrow(tab)>0) drawTable(tab[1:min(10,nrow(tab)),],x0=0,y0=0.98,dx=c(0.8,0.2),dy=0.04,row.names=FALSE,cex=0.5,col="#880000")
			
			par(fig=c(0.2,0.6,0.397,0.697),new=TRUE,mar=c(2,2,2,0));plot.new()
			tab = GO$GOMF[[icomp]]$neg[,c("Term","FDR")]
			tab = tab[tab$FDR<0.1,]
			tab$FDR =sprintf("%.2e",tab$FDR)
			text(0,1,sprintf("GO:MF neg : %d terms(FDR<0.01)",nrow(tab)),col="#000088",font=2,adj=c(0,0),cex=0.6)
			if (nrow(tab)>0) drawTable(tab[1:min(10,nrow(tab)),],x0=0,y0=0.98,dx=c(0.8,0.2),dy=0.04,row.names=FALSE,cex=0.5,col="#000088")
			par(fig=c(0.6,1,0.397,0.697),new=TRUE,mar=c(2,2,2,0));plot.new()
			tab = GO$GOMF[[icomp]]$pos[,c("Term","FDR")]
			tab = tab[tab$FDR<0.1,]
			tab$FDR =sprintf("%.2e",tab$FDR)
			text(0,1,sprintf("GO:MF pos : %d terms(FDR<0.01)",nrow(tab)),col="#880000",font=2,adj=c(0,0),cex=0.6)
			if (nrow(tab)>0) drawTable(tab[1:min(10,nrow(tab)),],x0=0,y0=0.98,dx=c(0.8,0.2),dy=0.04,row.names=FALSE,cex=0.5,col="#880000")
		}
		
		ix = 1; iy = 1
		if (!is.null(surv)){
			print("Survival")
			require("survival")
			par(fig=c(0.2+(ix-1)*0.2,0.2+ix*0.2,0.55-iy*0.2,0.55-(iy-1)*0.2),new=TRUE,mar=c(4,2,2,1))
			scoreD = c("low","high")[as.integer(IC$M[icomp,]>median(IC$M[icomp,]))+1]
			scoreD = factor(scoreD,levels=c("low","high"))
			cox = coxph(Surv(time = surv$time, event = surv$event) ~ IC$M[icomp,])
			pv = summary(cox)$logtest["pvalue"]
			lhr = log((summary(cox)$conf.int))[c(1,3,4)]
			if (pv<1e-2) {col=c("blue","red")}else{col=c("#666688","#886666")}
			plot(survfit(Surv(time = surv$time, event = surv$event) ~ scoreD,type="kaplan-meier"),col=col,conf.int=F,las=2,lwd=c(2,2,2),cex.axis=0.4)
			title(sprintf("Cox regression:\nlogtest pv=%.1e\nLHR=%.2f (CI = %.2f, %.2f)",pv,lhr[1],lhr[2], lhr[3]),cex.main=0.5)
		}
		
		if (!is.null(Var)){
			## ToDo: check whetehr factor is adequate (number of levels and NAs)
			print("ANOVA")
			pv = double(ncol(Var))+1
			names(pv)= names(Var)
			for (ifact in 1:ncol(Var)){
				fact = Var[[ifact]]
				ikeep = !is.na(fact)
				fact = as.factor(as.character(fact[ikeep]))
				if (nlevels(fact)==1) next
				x = IC$M[icomp,ikeep]
				res = aov(x~fact)
				pv[ifact] = summary(res)[[1]][1,5]
			}
			pv = pv[getTopIdx(-log10(pv),min(10,length(pv)))]
			tab = data.frame(factor = names(pv), p.value = sprintf("%.2e",pv))
			par(fig=c(0.4,1,0.4,0.58),new=TRUE,mar=c(0,2,2,1))
			plot.new()
			drawTable(tab,dx=c(0.8,0.2),dy=0.08,cex=0.5,row.names=FALSE,bg="white")
			ix=1; iy=2
			ifact=1
			for (ifact in 1:min(8,length(pv))){
				factname = names(pv)[ifact]
				fact = Var[[factname]]
				ikeep = !is.na(fact)
				fact = as.factor(as.character(fact[ikeep]))
				
				fontsize = 0.4
				if (nlevels(fact)<=5) fontsize = 0.6
				if (nlevels(fact)>15) fontsize = 0.3
				if (nlevels(fact)>30) fontsize = 0.2
				
				x = IC$M[icomp,ikeep]
				xf = list()
				for (i in 1:nlevels(fact)) xf[[i]] = x[fact == levels(fact)[i]]
				names(xf) = levels(fact)
				if (ix>4) {ix=1;iy=iy+1}
				par(fig=abs(c(0.2+(ix-1)*0.2,0.2+ix*0.2,0.6-iy*0.2,0.6-(iy-1)*0.2)),new=TRUE,mar=c(4,2,2,1))
				col="grey"
				if (pv[ifact]<1e-5)col="#66AAFF"
				if (pv[ifact]<1e-10)col="#88FF66"
				if (pv[ifact]<1e-20)col="#FF8866"
				violinplot(xf, col=col,cex.axis = fontsize,colMed = "black")
				title(sprintf("%s\npv=%.1e",factname,pv[ifact]),cex.main=0.5)
				ix = ix + 1
			}
		}
	}
	pheatmap(cor(t(IC$M))^2, main="R2 of M-matrix")
	if (!is.null(IC$stab)){
		par(mfcol=c(2,1),mar=c(4,3,3,1))
		boxplot(IC$stab,las=2,col="skyblue",main=sprintf("Stability (%d runs)",nrow(IC$stab)),subj="R2 between most correlates S in multiple runs")
	}
	if (!is.null(IC$mr2)) 
		plot(density(IC$mr2),lwd=2,col="blue",main="Distribution of mean R2\namong rows of M-matrix",xlab="Mean R2 for each single run")
	dev.off()
}

#!! ToDo: check whether mr2 proprly calculated!!


##=============================================================================
## Estimate optimal components for classification
## X - data matrix
## Y - vector of labels
##-----------------------------------------------------------------------------
estimateNComp = function (X,Y,kmin=2,kmax=ncol(X)%/%4,ntry=4,ncores=4){
	if (class(Y)!="factor") Y = factor(Y)
	cat(sprintf("Detection optimal numComponents for: %d features, %d samples (%d missing labels), %d classes, kmin=%d, kmax=%d\n",nrow(X),ncol(X),sum(is.na(Y)),nlevels(Y), kmin,kmax))
	Error = double(3)
	names(Error) = c("kmin","k","kmax")
	ic = runICA(X,ncomp=kmin,ntry=ntry,ncores=4)
	Error["kmin"] = runCrossValidation(t(ic$M),Y)$error
	ic = runICA(X,ncomp=kmax,ntry=ntry,ncores=4)
	Error["kmax"] = runCrossValidation(t(ic$M),Y)$error
	print(Error[-2])
	while (Error["kmax"] < Error["kmin"]){
		k = (kmax+kmin)%/%2
		ic = runICA(X,ncomp=k,ntry=ntry,ncores=4)
		Error["k"] = runCrossValidation(t(ic$M),Y)$error
		## worst case - E(k)>E(kmax) - shift k towards kmax
		if (Error["k"] > Error["kmax"]) {
			k = (kmax+k)%/%2
			Error["k"] = runCrossValidation(t(ic$M),Y)$error
		}
		cat(sprintf("k=%d, E(k)=%g, looking in [%d,%d] with E=(%g, %g)\n",k,Error["k"],kmin,kmax,Error["kmin"],Error["kmax"]))
		## shift borders
		if (Error["kmin"] > Error["k"]) {kmin = k; Error["kmin"]=Error["k"]}
		if (Error["kmax"] > Error["k"]) {kmax = k; Error["kmax"]=Error["k"]}
		if (kmax == kmin) {break}
	}
	return(k)
}

#estimateNComp(X, Var$Subtype) ## k is 60, looking in [2,119], error(k) = 0.0818731
#estimateNComp(X, Var$Sample.type) ## k is 60, looking in [2,119], error(k) = 0.127331

###############################################################################
## other functions used above
###############################################################################

##=============================================================================
## Sort dataframe, adapted from http://snippets.dzone.com/user/r-fanatic
##-----------------------------------------------------------------------------
sortDataFrame <- function(x, key, ...) {
    if (missing(key)) {
        rn <- rownames(x)
        if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
        x[order(rn, ...), , drop=FALSE]
    } else {
        x[do.call("order", c(x[key], ...)), , drop=FALSE]
    }
}
##=============================================================================
## num2fact - transforms numeric data to factors
##-----------------------------------------------------------------------------
num2fact = function(x, nlev = 4,digits=1){
	if (!class(x) %in% c("numeric","double","integer")) return(factor(x))
	f = rep(ifelse(digits==1,sprintf("q%d",1),sprintf("q%02d",1)),length(x))
	thr = quantile(x,seq(1/nlev,1,by = 1/nlev),na.rm=TRUE)
	for (ilev in 2:nlev)
		f[x > thr[ilev-1]] = ifelse(digits==1,sprintf("q%d",ilev),sprintf("q%02d",ilev))
	f[is.na(x)]=NA
	f=factor(f)
	return(f)
}
##=============================================================================
## df2fact - transforms data.frame to factors. Excludes factors with too small or large number of treatments
##-----------------------------------------------------------------------------
df2factor = function(Var,skip=TRUE,minlev =2, maxlev=Inf, NAs = c("","NA"), nlev=4, max.na.prop=0.9){
	if (class(Var)=="matrix") Var=as.data.frame(Var,stringsAsFactors=FALSE)
	if (!class(Var)=="data.frame") return(NULL)
	ikeep = logical(ncol(Var))|TRUE
	for (i in 1:ncol(Var)){
		if (sum(is.na(Var[[i]]))/nrow(Var) > max.na.prop) ikeep[i]= FALSE
		if (class(Var[[i]])%in%c("integer","numeric")) {
			Var[[i]] = num2fact(Var[[i]],nlev=nlev)
		}else{
			Var[[i]][Var[[i]] %in% NAs ] =NA
			Var[[i]] = factor(Var[[i]])
		}
		if (nlevels(Var[[i]])<minlev) ikeep[i]= FALSE
		if (nlevels(Var[[i]])>maxlev) ikeep[i]= FALSE
	}
	if (skip) Var=Var[,ikeep] 
	return(Var)
}

# ## Calculate confusion matrix
# getConfusionMatrix = function(gr, gr.pred) {
	# if (class(gr) != "factor") {gr = factor(gr); gr.pred = factor(gr.pred)}
	# nm = unique(levels(gr),levels(gr.pred))
	# Tab = matrix(nc = length(nm), nr = length(nm))
	# rownames(Tab) = paste("pred", nm, sep = ".")
	# colnames(Tab) = nm
	# for (i in 1:length(nm) )
		# for (j in 1:length(nm))
			# Tab[i,j] = sum((gr.pred == nm[i]) & (gr== nm[j]))
	# return(Tab)
# }

# # Computes the misclassification error from a confusion matrix.
# getMCError = function(CM) {
	# 1-sum(diag(CM))/sum(CM)
# }

# ## Computes accuracy
# getAccuracy = function(CM) {
	# sum(diag(CM))/sum(CM)
# }



################################################################################
## Tests
################################################################################
if (FALSE){
	setwd("b:/data/")
	download.file("http://data.sablab.net/astro.RData",mode="wb", destfile="astro.RData")
	load("astro.RData")
	str(X)
	X[,] = scale(X)
	
	groups = t(as.data.frame(strsplit(sub("_",".",colnames(X)),split="[.]")))
	rownames(groups) = colnames(X)
	groups = cbind(groups,sub("^.","",groups[,1]))
	groups = cbind(groups,sub(".$","",groups[,1]))
	colnames(groups) = c("Group","Patient","Platform","State","Tissue")

	ncomp = 8
	for (ncomp in c(2,4,6,8,12,16,24,32)){
		IC = runICA(X,ncomp=ncomp, ntry = 100, show.every=1, alg.typ="deflation",
					filter.thr = NULL, ncores=4, reduced=FALSE, savecomp=NULL, exclude=FALSE)
		GO = getGO(IC,alpha = 0.01)
		saveReport(IC,GO=GO, Var=groups,file = sprintf("report_ICA%d_noexcl.pdf",ncol(IC$S)))
		save(X,IC,groups,GO,file=sprintf("resICA%d_noexcl.RData",ncomp))
	}
	
	for (ncomp in c(2,4,6,8,12,16,24,32)){
		IC = runICA(X,ncomp=ncomp, ntry = 100, show.every=1, 
					filter.thr = NULL, ncores=4, reduced=FALSE, savecomp=NULL, exclude=TRUE)
		GO = getGO(IC,alpha = 0.01)
		saveReport(IC,GO=GO, Var=groups,file = sprintf("report_ICA%d_excl1.pdf",ncol(IC$S)))
		save(X,IC,groups,GO,file=sprintf("resICA%d_excl1.RData",ncomp))
	}
	
	write.table(IC$X,file="X.tsv",sep="\t",col.names=NA,quote=FALSE)
	write.table(IC$S,file="S.tsv",sep="\t",col.names=NA,quote=FALSE)
	write.table(IC$M,file="M.tsv",sep="\t",col.names=NA,quote=FALSE)
	write.table(groups,file="samples.tsv",sep="\t",col.names=NA,quote=FALSE)
	
	#ncomp=16;ntry = 10;show.every=1;filter.thr = 0;ncores=1;reduced=FALSE;savecomp=NULL;exclude=FALSE
	
	
	
	
	IC = runICA(X,ncomp=ncomp, ntry = 100, show.every=1, alg.typ="parallel",
				filter.thr = NULL, ncores=4, reduced=FALSE, savecomp=NULL, exclude=FALSE)
	GO = getGO(IC,alpha = 0.01)
	saveReport(IC,GO=GO, Var=groups,file = sprintf("report_ICA%d_noexcl_parallel.pdf",ncol(IC$S)))
	save(X,IC,groups,GO,file=sprintf("resICA%d_noexcl_parallel.RData",ncomp))
}