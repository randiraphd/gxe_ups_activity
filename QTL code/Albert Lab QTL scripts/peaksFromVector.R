# from a vector, make a list of indices that are higher than some threshold
# start a new vector whenever the entries are not continous
vectorAboveThreshold <- function(x, thres){
	ret = list(c())
	thresIndeces = which(x > thres)
# if nothing is above threshold, return empty list
	if (length(thresIndeces) == 0){
		return(ret)
	}
	if (length(thresIndeces) == 1){
		ret[[1]] = thresIndeces
		return(ret)
	}
	
	stepsizes = diff(thresIndeces)
	retFillIndex = 1
	for (i in 1:(length(thresIndeces)-1)){
		ret[[retFillIndex]] = c(ret[[retFillIndex]], thresIndeces[i])
		if (stepsizes[i] > 1){
			retFillIndex = retFillIndex + 1
			ret[[retFillIndex]] = as.integer(c())
		}
	}
# add the last element (retFillIndex will now automatically point to the correct index)
	ret[[retFillIndex]] = c(ret[[retFillIndex]], thresIndeces[length(thresIndeces)])
	ret
}

#test = vectorAboveThreshold(multiLOD[[12]][,2], 5)

# for a vector, find the positions of the local maxima:
# from http://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
# note that currently, this function also returns steps (e.g. 2 ,3 ,3, 4) => the first 3 and the 4 would be returned
# this is not ideal, although it might be that we can later fix it by combining the 2-LOD intervals?
localMaxima <- function(x) {
# Use -Inf instead if x is numeric (non-integer)
    if (length(unique(x)) == 1){return(1)}
	y <- diff(c(-Inf, x)) > 0
#	y <- diff(c(-.Machine$integer.max, x)) > 0L
	rle(y)$lengths
	y <- cumsum(rle(y)$lengths)
	y <- y[seq.int(1L, length(y), 2L)]
# this just seems to catch if the two first elements are the same?
    if (length(x) > 1){
        if (x[[1]] == x[[2]] & length(y) > 1) {
            y <- y[-1]
        }
    }
	y
}


# for a numeric vector and a vector of maxima (indeces), return a list of [somesize] unit drop intervals
# note that when using with maxima from a subregion of the whole chromosome, we need to add the start of that subregion!
xDropInterval <- function(LODs, maxima, dropsize){
	ret = matrix(ncol = 4, nrow = length(maxima))
	colnames(ret) = c("maxIndex", "maxValue", "leftIndex", "rightIndex")
	for (i in 1:length(maxima)){
		
		thisMax = LODs[maxima[i]]
		ret[i,] = c(maxima[i], thisMax, rep(maxima[i], 2))
		
# step left
		# currentPos will step away from the maximum out of 'maxima' that we currently look at
		currentPos = maxima[i]
		while(LODs[currentPos] >= (thisMax - dropsize) & currentPos > 1){
# make the step
			currentPos = currentPos - 1
		}
# once dropsize or the end of the chromomsome have been reached
		ret[i, 3] = currentPos

# step right
# currentPos will step away from the maximum out of 'maxima' that we currently look at
		currentPos = maxima[i]
		while(LODs[currentPos] >= (thisMax - dropsize) & currentPos < length(LODs)){
# make the step
			currentPos = currentPos + 1
		}
# once dropsize or the end of the chromomsome have been reached
		ret[i, 4] = currentPos
	}
	ret
}


# for a matrix of intervals, return, for each set of overlapping intervals, the one with the highest LOD
# (do not want to combine/extend intervals with overlapping LODdrops because of the localMaxima problem with steps)
# !!! assume that the matrix is in the format as produced by xDropInterval !!!
combineIntervals <- function(intervals){

	if (nrow(intervals) == 1){return(intervals)}

# make sure they are sorted	
	intervalsIntern = intervals[order(intervals[, 3]),]

	tempInterval = intervalsIntern[1,]
	ret = matrix(nrow = 1, ncol=4)
	retLine = 1
	ret[retLine,] = tempInterval
	for (i in 2:nrow(intervalsIntern)){
# if the intervals overlap, pick the current interval or tempInterval depending on who has the higher LOD
		if (tempInterval[4] >= intervalsIntern[i, 3]){
			tempInterval[3] = c(tempInterval[3], intervalsIntern[i, 3])[which.max(c(tempInterval[2], intervalsIntern[i, 2]))]
			tempInterval[4] = c(tempInterval[4], intervalsIntern[i, 4])[which.max(c(tempInterval[2], intervalsIntern[i, 2]))]
			tempInterval[1] = c(tempInterval[1], intervalsIntern[i, 1])[which.max(c(tempInterval[2], intervalsIntern[i, 2]))]
			tempInterval[2] = max(c(tempInterval[2], intervalsIntern[i, 2]))
			ret[retLine,] = tempInterval
# if they dont overlap, start a new resultinterval:
		}else{
			tempInterval = intervalsIntern[i, ]
			ret = rbind(ret, intervalsIntern[i, ])
			retLine = retLine + 1
		}
		
	}
    colnames(ret) = c("maxIndex", "maxValue", "leftIndex", "rightIndex")
	ret
}


#test2 = xDropInterval(multiLOD[[12]][,2], localMaxima(multiLOD[[12]][,2][test[[1]]]) + test[[1]][1], 2)
#xDropInterval(multiLOD[[12]][,2], which.max(multiLOD[[12]][,2][test[[1]]]) + test[[1]][1], 2)

#combineIntervals(test2)
#test3 = test2
#test3[,1] = test3[,1] + 10e3
#test3[,2] = test3[,2] + 2
#test3[,3] = test3[,3] + 10e3
#test3[,4] = test3[,4] + 10e3
#test4 = rbind(test2, test3)
#combineIntervals(test4)


callPeaks <- function(LOD, thres, dropSize){
    hiLODRegions = vectorAboveThreshold(LOD, thres)
    if (is.null(hiLODRegions[[1]])){return(c())}
    LODIntervals = c()
    for (x in hiLODRegions){
        LODIntervals = rbind(
            LODIntervals, 
            combineIntervals(xDropInterval(LOD, localMaxima(LOD[x]) + x[1] - 1, dropSize))
#			combineIntervals(xDropInterval(LOD, which(LOD[x] == max(LOD[x])) + x[1] - 1, dropSize))				 
        )
    }
    # this extra combination step is necessary because sometimes two close peaks can be split into two hiLODregions if the dip goes below threshold
    LODIntervals = combineIntervals(LODIntervals)
    # these adjustments bring the coordinates into bp space
    # note that 100bp was my binsize for running multipool
    LODIntervals[,1] = (LODIntervals[,1]-1) * 100
    LODIntervals[,3] = (LODIntervals[,3]-1) * 100
    LODIntervals[,4] = (LODIntervals[,4]-1) * 100
    LODIntervals
}

#callPeaks(multiLOD[[12]][,2], 5, 2)

