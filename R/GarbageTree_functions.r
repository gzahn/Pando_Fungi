# monster_tree_functions.r
require(data.table)
require(phytools)
require(ape)



# function to read in fasta or fastq file as list
# NOTE: ONLY WORKS ON SEQUENTIAL FASTX FILES. 
# Hope you aren't still using interleaved files after like 2005
read.fastx <- function(file, type=c("fastq", "fasta")){
	# read in raw data
	lines <- scan(file, what="character", sep='\n')
	
	# check to make sure data are OK-ish
	nlines <- length(lines)
	if(type=="fastq" & nlines %% 4 > 0){
		stop("CRITICAL ERROR: number of lines in fastq file not divisible by 4.")
	}else if(type == "fasta" & nlines %% 2 > 0){
		stop("CRITICAL ERROR: number of lines in fasta file not divisible by 2")
	}

	# make indices for item types (1=header, 2=seq, 3="+", 4=qual)
	if(type=="fastq"){
		line_inds <- rep(1:4, nlines/4)
	}else if(type=="fasta"){
		line_inds <- rep(1:2, nlines/2)
	}

	# make names generic
	headers <- as.vector(sapply(X=lines[line_inds == 1], FUN=function(x) substring(x, 2) ) )

	# make output list object
	if(type=="fastq"){
		output <- mapply(FUN=c, lines[line_inds == 2], lines[line_inds == 4], SIMPLIFY=FALSE)
	}else if(type=="fasta"){
		output <- mapply(FUN=c, lines[line_inds == 2], SIMPLIFY=FALSE)
	}

	# add names to output list
	names(output) <- headers

	# all done
	return(output)

}

# turns the above list into a simple string array of headers and seqs
# that can be written to text as a fasta file
fastalist2char <- function(fastalist){
	outvec <- rep("error", length(fastalist) * 2)
	for(i in 1:length(fastalist)){
		writepos <- (i * 2) - 1
		# write header
		outvec[writepos] <- paste(">", names(fastalist)[i], sep="")
		# write seq
		outvec[writepos+1] <- fastalist[[i]][1]
	}
	return(outvec)
}


# runs vsearch on sequences
vsearch_r <- function(seqs, id_cutoff, tmpfile_prefix="vsearch_tmp_seqs"){
	tmpfile <- paste(tmpfile_prefix, "_", Sys.getpid(), ".fasta", sep="")
	require(data.table)
	# write seqs to tmpfile
	# the "input" argument for system() is NOT working properly...
	fwrite(list(fastalist2char(seqs)), file=tmpfile)

	# make vsearch string
	v_string <- paste("vsearch --cluster_fast ", tmpfile, " --id", id_cutoff, "--uc -")

	# run vsearch, capture uc output!
	uctable <- system(command=v_string, intern=T, ignore.stderr=T)
	uctable <- read.table(text=uctable, sep='\t', header=F, stringsAsFactors=F)

	# remove temp file
	rmstring <- paste("rm", tmpfile)
	system(command=rmstring)

	return(uctable)
}

# runs vsearch iteratively at different clustering levels.
# returns list of uc tables
iterative_cluster <- function(otuseqs, startlvl=0.99, finallvl=0.10, interval=0.02){
	# a vector of cutoffs to iterate over
	cutoffs <- seq(from=startlvl, to=finallvl, by=-interval)
	# initialize list object
	clustering_tables_list <- list()
	# set up progress bar
	pb <- txtProgressBar(min=0, max=length(cutoffs), style=3)
	# iteratively cluster
	for(i in 1:length(cutoffs)){
		# subset sequences and rename if necessary
		if(i == 1){
			seqs_i <- otuseqs
		}else{
			seqs2keep <- unique(clustering_tables_list[[i-1]]$new)
			seqs_i <- otuseqs[names(otuseqs) %in% seqs2keep]
		}
		# run vsearch
		uctable_i <- vsearch_r(seqs=seqs_i, id_cutoff=cutoffs[i])
		# prune useless centroid records
		uctable_i <- uctable_i[uctable_i[[1]] != "C", ]
		# make centroids reference themselves
		c_rows <- (uctable_i[[10]] == "*")
		uctable_i[[10]][c_rows] <- uctable_i[[9]][c_rows]
		# dump useless rows in uc table
		uctable_i <- data.frame(old=uctable_i[[9]], new=uctable_i[[10]])
		# add result to clustering_tables_list
		clustering_tables_list[[i]] <- uctable_i
		names(clustering_tables_list)[[i]] <- paste("id=", cutoffs[i], sep="")
		# check if everything is 100% clustered. if so, break
		if(all(uctable_i$new == uctable_i$new[1])){
			setTxtProgressBar(pb, length(cutoffs))
			print("Reached max clustering early.")
			break()
		}
		# update progress bar
		setTxtProgressBar(pb, i)
	}

	return(clustering_tables_list)


}

# turns each entry in a uc table into a nwk string for a tip
# returns a table of nwk strings for each tip
nwkize_tips <- function(uc, delta){
	centroids <- unique(uc$new)
	nwk_strs <- rep("", length(centroids))
	for(i in 1:length(centroids)){
		tips_cent_i <- uc$old[ uc$new == centroids[i] ]
		tips_cent_i <- paste(tips_cent_i, ":", delta, sep="")
		tips_cent_i <- paste(tips_cent_i, collapse=",")
		nwk_strs[i] <- paste("(", tips_cent_i, ")", sep="")
	}
	return(data.frame(centroids, nwk_strs, stringsAsFactors=FALSE))
}

# combines nwk strings together per a uc table
# takes output of this function or nwkize_tips as input
# returns a table of nwk strings for each node
nwkize_nodes <- function(nwks_tab, uc, delta){
	centroids <- unique(uc$new)
	nwk_strs <- rep("", length(centroids))
	for(i in 1:length(centroids)){
		tips_cent_i <- uc$old[ uc$new == centroids[i] ]
		nwks_cent_i <- nwks_tab$nwk_strs[ nwks_tab$centroids %in% tips_cent_i ]
		#if(length(nwks_cent_i) == 1){
		#	colon_positions <- as.numeric( gregexpr(pattern =':',nwks_cent_i[1])[[1]] )
		#	last_colon_pos <- colon_positions[length(colon_positions)]
		#	before <- substr(nwks_cent_i[1], start=2, stop=(last_colon_pos - 1))
		#	last_delta <- as.numeric(substr(nwks_cent_i[1], start=(last_colon_pos + 1), stop=(nchar(nwks_cent_i[1]) - 1)))
		#	nwk_strs[i] <- paste("(", before, ":", delta+last_delta, ")", sep="")
		#}else{
			nwks_cent_i <- paste(nwks_cent_i, ":", delta, sep="")
			nwks_cent_i <- paste(nwks_cent_i, collapse=",")
			nwk_strs[i] <- paste("(", nwks_cent_i, ")", sep="")
		#}
	}
	return(data.frame(centroids, nwk_strs, stringsAsFactors=FALSE))
}

# wrapper script for nwkize_tips() and nwkize_nodes().
# turns a list of uc tables into a cladogram.
# note: uc table list must be sorted, from tips to root.
# this is the default output of iterative_cluster()
make_garbage_tree <- function(uc_list, interval=-1, deltas=-1){
	# note - this requires tips to have low indices, and basal nodes to have high indices
	# ALL OBECTS MUST BE SORTED HI->LO in terms of TIPS->ROOT
	# build deltas unless custom deltas are specified
	if(interval == -1 && deltas[1] == -1){
		stop("ERROR: no interval specified")
	}else if(interval > 0){
		deltas <- rep(interval, length(uc_list))
	}
	# progress message:
	message("Building tree from UC tables:")
	# progress is proportional to the number of rows to collapse:
	total_rows <- sum(unlist(lapply(FUN=nrow, X=uc_list)))
	rows_done <- 0
	pb <- txtProgressBar(min=0, max=total_rows, style=3)
	# make tree, top-down
	for(i in 1:length(uc_list)){
		if(i == 1){
			last_nwks_tab <- nwkize_tips(uc=uc_list[[i]], delta=deltas[i])
		}else{
			last_nwks_tab <- nwkize_nodes(nwks_tab=last_nwks_tab, uc=uc_list[[i]], delta=deltas[i])		
		}
		rows_done <- rows_done + nrow(uc_list[[i]])
		setTxtProgressBar(pb, rows_done)
	}
	# create basal polytomy
	if(nrow(last_nwks_tab) > 1){
		tree <- paste(last_nwks_tab$nwk_strs, collapse=":0,")
		tree <- paste("(", tree, ":0)", sep="")
	}else{
		tree <- last_nwks_tab$nwk_strs[1]
	}
	tree <- paste(tree, ":0;", sep="")
	message("...done.")
	# convert tree to R phylo object
	tree <- read.tree(text=tree)

	tree <- collapse.singles(tree)


	# fix NA branch lengths (should be zero)
	# tree$edge.length[is.na(tree$edge.length)] <- 0

	return(tree)
}

simplify_garbage_nodes <- function(tree){
	return(collapse.singles(tree))
}
 