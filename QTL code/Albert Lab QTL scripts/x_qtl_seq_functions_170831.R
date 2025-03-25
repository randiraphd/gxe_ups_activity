# this version attempts at dynamic smoothing for different chr lengths (otherwise short chromosomes get very jumpy)
# in practice, we then switched to using coverage as the weight???
rollLoessByChrWithWeights = function(x, span4Loess=0.1){
	unlist(sapply(unique(x[,1]), function(y){
				  temp = x[x[,1] == y,]
				  ret = rep(NA, nrow(temp))
#				  print(c(length(temp[,3][is.na(temp[,3])]), nrow(temp)))
				  if(length(temp[,3][is.na(temp[,3])]) < nrow(temp)){
				  try({
					  ret = predict(loess(as.numeric(temp[,2]) ~ as.numeric(temp[,3]), weights = as.numeric(temp[,4]), span=1/((span4Loess/40) * nrow(temp[complete.cases(temp),]))), as.numeric(temp[,3]))
					  })
				  }
				  names(ret) = rownames(temp)
				  ret
				  }))
}

getGcoords = function ( chr , pos, spacing=0, sgd.table=paste0(proj_dir, "QTL_scripts/sacCer3ChromLengths.txt" )) {
    offind = as.vector(cumsum(read.table(sgd.table, header=FALSE, sep="\t")[,2] + spacing))
    offind=    offind[-length(offind)]
    offind= c(0, offind)
    names(offind) = as.character(read.table(sgd.table, header=FALSE, sep="\t")[,1])
    chr.off=as.numeric(sapply(chr, function(x) {offind[[x]]}))
    return(chr.off+pos)
}
