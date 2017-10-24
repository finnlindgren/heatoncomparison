#
# Source file for agglomerative and divisive clustering algorithm
# Documentation in provided
#

library(LatticeKrig)
library(deldir)

div.gradclust <- function(y,s,num.clust,bndry.pts,method="ward",num.gpts=NULL){
	
	if(num.clust<=0){
		
		stop("Inappropriate number of clusters: 1<= num.clust <= nrow(s)")
		
	} else if(num.clust>nrow(s)) {
		
		stop("Inappropriate number of clusters: 1<= num.clust <= nrow(s)")
		
	} else if(num.clust==1){
		
		## If num.clust==1 return only 1 cluster
		clust <- rep(1,nrow(s))
		
	} else {
		
		if(is.null(num.gpts)){
			
			## if num.gpts is null then cluster the observations
			cat("num.gpts not provided...Clustering raw observations...\n")
			clust <- obs.div.gradclust(y=y, s=s, num.clust=num.clust, bndry.pts=bndry.pts,method=method)
			
		} else {
			
			## else run a divisive clustering on an aggregated lattice
			cat("num.gpts provided...Clustering on lattice... \n")
			clust <- grid.div.gradclust(y=y, s=s, num.clust=num.clust, bndry.pts=bndry.pts, method=method, num.gpts=num.gpts)
			
		}
		
	}
	
	return(clust)
	
}

obs.div.gradclust <- function(y,s,num.clust,bndry.pts,method){
	
	## Get the Voronoi Tesselation
	cat("Determining Voronoi Tesselation of s...\n")
	v.tess <- deldir(s[,1],s[,2],dpl=list(ndx=2,ndy=2),rw=bndry.pts)
	neighbors <- vector("list",nrow(s))
	distances <- vector("list",nrow(s))
	changes <- vector("list",nrow(s))
	for(loc in 1:nrow(s)){
		the.rows <- which(((v.tess$dirsgs[,5]==loc)+(v.tess$dirsgs[,6]==loc))>0)
		the.neighbors <- c(as.matrix(v.tess$dirsgs[the.rows,5:6]))
	
		neighbors[[loc]] <- the.neighbors[the.neighbors!=loc]
		neighbors[[loc]] <- neighbors[[loc]][neighbors[[loc]]<=nrow(s)]
		
		distances[[loc]] <- as.numeric(rdist(matrix(s[loc,],nrow=1),matrix(s[neighbors[[loc]],],ncol=2)))
		#changes[[loc]] <- (abs(rep(y[loc],length(neighbors[[loc]]))-y[neighbors[[loc]]]))
		changes[[loc]] <- (abs(rep(y[loc],length(neighbors[[loc]]))-y[neighbors[[loc]]]))/(distances[[loc]])
	}
	
	## Do Divisive Clustering
	cat("Determining Cluster Membership - can be time consuming for large nrow(s)...\n")
	clust <- rep(1,nrow(s))
	for(k in 1:(num.clust-1)){
		
		## Find Largest Change within each cluster
		possible.breaks <- c()
		for(i in sort(unique(clust))){
			#Determine Points in Cluster
			clust.pts <- which(clust==i)
			
			## Go through each point in the cluster and find maximum gradient within cluster
			if(length(clust.pts)>1){
				max.break <- c() 
				for(j in clust.pts){
					## Within cluster Neighbors/Changes of j
					chg <- changes[[j]]
					nn <- neighbors[[j]]
					ds <- distances[[j]]
					chg <- chg[clust[nn]==i]
					ds <- ds[clust[nn]==i]
					nn <- nn[clust[nn]==i]
					
					## Keep Largest within cluster match
					max.break <- rbind(max.break,cbind(j,nn[chg==max(chg)],chg[chg==max(chg)]))
					
				}
			} else {
				## Don't split a cluster with only 1 obs in it
				max.break <- matrix(c(clust.pts,clust.pts+1,-1),nrow=1)
			}
		
			## Keep only maximum break within each cluster
			possible.breaks <- rbind(possible.breaks,max.break[max.break[,3]==max(max.break[,3]),])
		
		} #End i loop
		possible.breaks <- matrix(possible.breaks[possible.breaks[,1]<possible.breaks[,2],],ncol=3)
	
		## Select break.pt and reassign label
		split.clust <- which(possible.breaks[,3]==max(possible.breaks[,3]))
		tmp.assigned <- rep(0,nrow(s))
		tmp.assigned[clust!=split.clust] <- 1
		break.pt <- possible.breaks[split.clust,1:2]
		clust[break.pt[2]] <- k+1
	
		## Reassign remaining points to a cluster from closest to break to farthest away
		tmp.assigned[break.pt] <- 1
		while(sum(tmp.assigned)<nrow(s)){
			## Get Neighboring Points of Previously Assigned
			neigh <- c()
			for(i in which((tmp.assigned==1)+(clust%in%c(split.clust,k+1))==2)){
				neigh <- rbind(neigh,cbind(rep(clust[i],length(neighbors[[i]])),neighbors[[i]]))
			}
			neigh <- matrix(neigh[which((tmp.assigned[neigh[,2]]==0)+(clust[neigh[,2]]%in%c(split.clust,k+1))==2),],ncol=2)
		
			## Just Assign Points that Neighbor >1 Clusters
			num.neigh <- table(neigh[,2],neigh[,1])
			neigh <- row.names(num.neigh)
			neigh <- as.numeric(neigh[rowSums(num.neigh>0)==max(rowSums(num.neigh>0))])
					
			## Do the Assignment
			cmns <- aggregate(y[tmp.assigned==1],by=list(clust=clust[tmp.assigned==1]),FUN=mean)
			clust.size <- aggregate(tmp.assigned[tmp.assigned==1],by=list(cluster=clust[tmp.assigned==1]),FUN=sum)
			for(i in neigh){
				nn <- neighbors[[i]]
				chg <- changes[[i]]
				ds <- distances[[i]]
				the.chgs <- abs(cmns[clust[nn],2]-y[i])/distances[[i]]
				chg <- chg[nn%in%which(tmp.assigned==1)]
				the.chgs <- the.chgs[nn%in%which(tmp.assigned==1)]
				ds <- ds[nn%in%which(tmp.assigned==1)]
				nn <- nn[nn%in%which(tmp.assigned==1)]
				chg <- chg[clust[nn]%in%c(split.clust,k+1)]
				the.chgs <- the.chgs[clust[nn]%in%c(split.clust,k+1)]
				ds <- ds[clust[nn]%in%c(split.clust,k+1)]
				nn <- nn[clust[nn]%in%c(split.clust,k+1)]
				if(method=="ward"){
					ward <- (clust.size[clust[nn],2]/(1+clust.size[clust[nn],2]))*((cmns[clust[nn],2]-y[i])^2)/ds
					clust[i] <- clust[nn[ward==min(ward)]]
				} else if(method=="complete"){
					clust[i] <- unique(clust[nn[the.chgs==min(the.chgs)]])
				} else if(method=="single"){
					clust[i] <- clust[nn[chg==min(chg)]]
				}
				tmp.assigned[i] <- 1
			}
			
		} #End while loop
		
	} # End k loop
	cat('Clustering complete.\n')
	return(clust)
} #End obs.div.gradclust

grid.div.gradclust <- function(y,s,num.clust,num.gpts,bndry.pts,method){
	
	## Define a Grid on which to do agglomerative clustering
	xlim <- range(s[,1])
	ylim <- range(s[,2])
	xseq <- seq(xlim[1],xlim[2],length=num.gpts+2)
	yseq <- seq(ylim[1],ylim[2],length=num.gpts+2)
	gpts <- as.matrix(expand.grid(xseq[2:(length(xseq)-1)],yseq[2:(length(yseq)-1)]))
	
	## Aggregate Observations up to grid points
	dists.gpts <- rdist(gpts,s)
	assigned.gpt <- apply(dists.gpts==matrix(apply(dists.gpts,2,min),nrow=nrow(dists.gpts),ncol=nrow(s),byrow=TRUE),2,which)
	mn.gpts <- aggregate(y,by=list(gpt=assigned.gpt),FUN=mean)
	gpts <- gpts[mn.gpts[,1],]
	mn.gpts <- mn.gpts[,2]
	dists.gpts <- rdist(gpts,s)
	assigned.gpt <- apply(dists.gpts==matrix(apply(dists.gpts,2,min),nrow=nrow(dists.gpts),ncol=nrow(s),byrow=TRUE),2,which)
	
	## Get the Voronoi Tesselation
	cat("Determining Voronoi Tesselation...\n")
	v.tess <- deldir(gpts[,1],gpts[,2],dpl=list(ndx=2,ndy=2),rw=bndry.pts)
	neighbors <- vector("list",nrow(gpts))
	distances <- vector("list",nrow(gpts))
	changes <- vector("list",nrow(gpts))
	for(loc in 1:nrow(gpts)){
		the.rows <- which(((v.tess$dirsgs[,5]==loc)+(v.tess$dirsgs[,6]==loc))>0)
		the.neighbors <- c(as.matrix(v.tess$dirsgs[the.rows,5:6]))
	
		neighbors[[loc]] <- the.neighbors[the.neighbors!=loc]
		neighbors[[loc]] <- neighbors[[loc]][neighbors[[loc]]<=nrow(gpts)]
		
		distances[[loc]] <- as.numeric(rdist(matrix(gpts[loc,],nrow=1),matrix(gpts[neighbors[[loc]],],ncol=2)))
		#changes[[loc]] <- (abs(rep(y[loc],length(neighbors[[loc]]))-y[neighbors[[loc]]]))
		changes[[loc]] <- (abs(rep(mn.gpts[loc],length(neighbors[[loc]]))-mn.gpts[neighbors[[loc]]]))/(distances[[loc]])
	}
	
	## Do Divisive Clustering
	cat("Determining Cluster Membership - can be time consuming for large num.gpts...\n")
	clust <- rep(1,nrow(gpts))
	for(k in 1:(num.clust-1)){
		
		## Find Largest Change within each cluster
		possible.breaks <- c()
		for(i in sort(unique(clust))){
			#Determine Points in Cluster
			clust.pts <- which(clust==i)
			
			## Go through each point in the cluster and find maximum gradient within cluster
			if(length(clust.pts)>1){
				max.break <- c() 
				for(j in clust.pts){
					## Within cluster Neighbors/Changes of j
					chg <- changes[[j]]
					nn <- neighbors[[j]]
					ds <- distances[[j]]
					chg <- chg[clust[nn]==i]
					ds <- ds[clust[nn]==i]
					nn <- nn[clust[nn]==i]
					
					## Keep Largest within cluster match
					max.break <- rbind(max.break,cbind(j,nn[chg==max(chg)],chg[chg==max(chg)]))
				}
			} else {
				## Don't split a cluster with only 1 obs in it
				max.break <- matrix(c(clust.pts,clust.pts+1,-1),nrow=1)
			}
		
			## Keep only maximum break within each cluster
			possible.breaks <- rbind(possible.breaks,max.break[max.break[,3]==max(max.break[,3]),])
		
		} #End i loop
		possible.breaks <- matrix(possible.breaks[possible.breaks[,1]<possible.breaks[,2],],ncol=3)
	
		## Select break.pt and reassign label
		split.clust <- which(possible.breaks[,3]==max(possible.breaks[,3]))
		tmp.assigned <- rep(0,nrow(gpts))
		tmp.assigned[clust!=split.clust] <- 1
		break.pt <- possible.breaks[split.clust,1:2]
		clust[break.pt[2]] <- k+1
	
		## Reassign remaining points to a cluster from closest to break to farthest away
		tmp.assigned[break.pt] <- 1
		while(sum(tmp.assigned)<nrow(gpts)){
			## Get Neighboring Points of Previously Assigned
			neigh <- c()
			for(i in which((tmp.assigned==1)+(clust%in%c(split.clust,k+1))==2)){
				neigh <- rbind(neigh,cbind(rep(clust[i],length(neighbors[[i]])),neighbors[[i]]))
			}
			neigh <- matrix(neigh[which((tmp.assigned[neigh[,2]]==0)+(clust[neigh[,2]]%in%c(split.clust,k+1))==2),],ncol=2)
		
			## Just Assign Points that Neighbor >1 Clusters
			num.neigh <- table(neigh[,2],neigh[,1])
			neigh <- row.names(num.neigh)
			neigh <- as.numeric(neigh[rowSums(num.neigh>0)==max(rowSums(num.neigh>0))])
					
			## Do the Assignment
			cmns <- aggregate(mn.gpts[tmp.assigned==1],by=list(clust=clust[tmp.assigned==1]),FUN=mean)
			clust.size <- aggregate(tmp.assigned[tmp.assigned==1],by=list(cluster=clust[tmp.assigned==1]),FUN=sum)
			for(i in neigh){
				nn <- neighbors[[i]]
				chg <- changes[[i]]
				ds <- distances[[i]]
				the.chgs <- abs(cmns[clust[nn],2]-y[i])/distances[[i]]
				chg <- chg[nn%in%which(tmp.assigned==1)]
				the.chgs <- the.chgs[nn%in%which(tmp.assigned==1)]
				ds <- ds[nn%in%which(tmp.assigned==1)]
				nn <- nn[nn%in%which(tmp.assigned==1)]
				chg <- chg[clust[nn]%in%c(split.clust,k+1)]
				the.chgs <- the.chgs[clust[nn]%in%c(split.clust,k+1)]
				ds <- ds[clust[nn]%in%c(split.clust,k+1)]
				nn <- nn[clust[nn]%in%c(split.clust,k+1)]
				if(method=="ward"){
					ward <- (clust.size[clust[nn],2]/(1+clust.size[clust[nn],2]))*((cmns[clust[nn],2]-mn.gpts[i])^2)/ds
					clust[i] <- unique(clust[nn[ward==min(ward)]])
				} else if(method=="complete"){
					clust[i] <- unique(clust[nn[the.chgs==min(the.chgs)]])
				} else if(method=="single"){
					clust[i] <- clust[nn[chg==min(chg)]]
				}
				tmp.assigned[i] <- 1
			}
			
		} #End while loop
		
	} # End k loop
	clust <- clust[assigned.gpt]
	cat('Clustering complete.\n')
	return(clust)
} # End grid.div.gradclust

agg.gradclust <- function(y,s,num.clust,bndry.pts,method="ward",num.gpts=NULL){
	
	if(num.clust<=0){
		
		stop("Inappropriate number of clusters: 1<= num.clust <= nrow(s)")
		
	} else if(num.clust>nrow(s)) {
		
		stop("Inappropriate number of clusters: 1<= num.clust <= nrow(s)")
		
	} else if(num.clust==1){
		
		## If num.clust==1 return only 1 cluster
		clust <- rep(1,nrow(s))
		
	} else {
		
		if(is.null(num.gpts)){
			
			## if num.gpts is null then cluster the observations
			cat("num.gpts not provided: clustering raw observations...\n")
			clust <- obs.agg.gradclust(y=y, s=s, num.clust=num.clust, bndry.pts=bndry.pts,method=method)
			
		} else {
			
			## else run a divisive clustering on an aggregated lattice
			cat("num.gpts provided: clustering on lattice...\n")
			clust <- grid.agg.gradclust(y=y, s=s, num.clust=num.clust, bndry.pts=bndry.pts, method=method, num.gpts=num.gpts)
			
		}
		
	}
	
	return(clust)
	
} ## End agg.gradclust

obs.agg.gradclust <- function(y,s,num.clust,bndry.pts,method="ward"){
	
	## Get the Voronoi Tesselation
	cat("Determining Voronoi Tesselation...\n")
	v.tess <- deldir(s[,1],s[,2],dpl=list(ndx=2,ndy=2),rw=bndry.pts)
	neighbors <- vector("list",nrow(s))
	distances <- vector("list",nrow(s))
	changes <- vector("list",nrow(s))
	all.changes <- c()
	for(loc in 1:nrow(s)){
		the.rows <- which(((v.tess$dirsgs[,5]==loc)+(v.tess$dirsgs[,6]==loc))>0)
		the.neighbors <- c(as.matrix(v.tess$dirsgs[the.rows,5:6]))
	
		neighbors[[loc]] <- the.neighbors[the.neighbors!=loc]
		neighbors[[loc]] <- neighbors[[loc]][neighbors[[loc]]<=nrow(s)]
		
		distances[[loc]] <- as.numeric(rdist(matrix(s[loc,],nrow=1),s[neighbors[[loc]],]))
		changes[[loc]] <- (abs(rep(y[loc],length(neighbors[[loc]]))-y[neighbors[[loc]]])^2)/distances[[loc]]
		all.changes <- rbind(all.changes,cbind(changes[[loc]],rep(loc,length(neighbors[[loc]])),neighbors[[loc]],distances[[loc]]))
	}
	all.changes <- all.changes[,c(2,3,4,1)]
	
	## Do Agglomerative Clustering
	cat("Determining Cluster Membership - can be time consuming for large nrow(s)...\n")
	clust <- 1:length(y)
	for(k in 1:(length(y)-num.clust)){
		
		## Clock Expected Computation time
		if(k==1){
			t0 <- proc.time()[1:3]
		}
		
		## Calculate Dissimilarity Between Each Cluster
		if(method!="ward"){
			possible.mergers <- cbind(clust[all.changes[,1]],clust[all.changes[,2]],all.changes[,4])
			diff.clust <- possible.mergers[,1]!=possible.mergers[,2]
			possible.mergers <- possible.mergers[diff.clust,]
			if(method=="complete"){
				possible.mergers <- aggregate(possible.mergers[,3],by=list(clstr1=possible.mergers[,1],clstr2=possible.mergers[,2]),FUN=max)
			} else if(method=="single"){
				possible.mergers <- aggregate(possible.mergers[,3],by=list(clstr1=possible.mergers[,1],clstr2=possible.mergers[,2]),FUN=min)
			} else if(method=="average"){
				possible.mergers <- aggregate(possible.mergers[,3],by=list(clstr1=possible.mergers[,1],clstr2=possible.mergers[,2]),FUN=mean)
			}
			possible.mergers <- possible.mergers[order(possible.mergers[,1]),]
			possible.mergers <- as.matrix(possible.mergers[possible.mergers[,1]<possible.mergers[,2],])
		
		} else {
		
			## Calculate Ward Dissimilarity Between Each Cluster
			cmns <- aggregate(y,by=list(cluster=clust),FUN=mean)
			possible.mergers <- aggregate(all.changes[,3],by=list(c1=clust[all.changes[,1]],c2=clust[all.changes[,2]]),FUN=mean)
			possible.mergers <- as.matrix(possible.mergers[possible.mergers[,1]<possible.mergers[,2],])
			clust.sizes <- table(clust)
			ward <- ((cmns[possible.mergers[,2],2]-cmns[possible.mergers[,1],2])^2)*(clust.sizes[possible.mergers[,1]]*clust.sizes[possible.mergers[,2]])/(clust.sizes[possible.mergers[,1]]+clust.sizes[possible.mergers[,2]])
			possible.mergers <- as.matrix(cbind(possible.mergers[,1:2],ward/possible.mergers[,3]))
		
		}
		
		## Determine best merger = smallest dissimilarity
		if(!is.null(dim(possible.mergers))){
			merge.pts <- possible.mergers[possible.mergers[,3]==min(possible.mergers[,3]),1:2]
		} else {
			merge.pts <- possible.mergers[possible.mergers[3]==min(possible.mergers[3])]
		}
	
		## Merge and Recalculate number of clusters
		new.lab <- min(c(merge.pts[1],merge.pts[2]))
		clust[clust==merge.pts[1]] <- new.lab
		clust[clust==merge.pts[2]] <- new.lab
		clust[clust>max(merge.pts)] <- clust[clust>max(merge.pts)]-1
		
		## Print out expected computation time
		if(k==1){
			t1 <- proc.time()[1:3]-t0
			cat('estimated clustering time (min):',t1[3]*(length(y)-num.clust)/60,'...\n')
		}
	} # End k loop
	clust <- as.numeric(factor(clust))
	cat('Clustering complete.\n')
	return(clust)
	
} #End obs.agg.gradclust

grid.agg.gradclust <- function(y,s,num.clust,num.gpts,bndry.pts,method="ward"){
	
	## Define a Grid on which to do agglomerative clustering
	xlim <- range(s[,1])
	ylim <- range(s[,2])
	xseq <- seq(xlim[1],xlim[2],length=num.gpts+2)
	yseq <- seq(ylim[1],ylim[2],length=num.gpts+2)
	gpts <- as.matrix(expand.grid(xseq[2:(length(xseq)-1)],yseq[2:(length(yseq)-1)]))
	
	## Aggregate Observations up to grid points
	dists.gpts <- rdist(gpts,s)
	assigned.gpt <- apply(dists.gpts==matrix(apply(dists.gpts,2,min),nrow=nrow(dists.gpts),ncol=nrow(s),byrow=TRUE),2,which)
  if(class(assigned.gpt)=="list"){
    assigned.gpt <- unlist(lapply(assigned.gpt,function(x){return(x[1])}))
  }
	mn.gpts <- aggregate(y,by=list(gpt=assigned.gpt),FUN=mean)
	gpts <- gpts[mn.gpts[,1],]
	mn.gpts <- mn.gpts[,2]
	dists.gpts <- rdist(gpts,s)
	assigned.gpt <- apply(dists.gpts==matrix(apply(dists.gpts,2,min),nrow=nrow(dists.gpts),ncol=nrow(s),byrow=TRUE),2,which)
	if(class(assigned.gpt)=="list"){
	  assigned.gpt <- unlist(lapply(assigned.gpt,function(x){return(x[1])}))
	}

	## Get the Voronoi Tesselation
	cat("Determining Voronoi Tesselation...\n")
	v.tess <- deldir(gpts[,1],gpts[,2],dpl=list(ndx=2,ndy=2),rw=bndry.pts)
	neighbors <- vector("list",nrow(gpts))
	distances <- vector("list",nrow(gpts))
	changes <- vector("list",nrow(gpts))
	all.changes <- c()
	for(loc in 1:nrow(gpts)){
		the.rows <- which(((v.tess$dirsgs[,5]==loc)+(v.tess$dirsgs[,6]==loc))>0)
		the.neighbors <- c(as.matrix(v.tess$dirsgs[the.rows,5:6]))
	
		neighbors[[loc]] <- the.neighbors[the.neighbors!=loc]
		neighbors[[loc]] <- neighbors[[loc]][neighbors[[loc]]<=nrow(gpts)]
		
		distances[[loc]] <- as.numeric(rdist(matrix(gpts[loc,],nrow=1),gpts[neighbors[[loc]],]))
		changes[[loc]] <- (abs(rep(mn.gpts[loc],length(neighbors[[loc]]))-mn.gpts[neighbors[[loc]]])^2)/distances[[loc]]
		all.changes <- rbind(all.changes,cbind(changes[[loc]],rep(loc,length(neighbors[[loc]])),neighbors[[loc]],distances[[loc]]))
	}
	all.changes <- all.changes[,c(2,3,4,1)]
	
	## Do Agglomerative Clustering
	cat("Determining Cluster Membership - can be time consuming for large num.gpts...\n")
	clust <- 1:length(mn.gpts)
	for(k in 1:(length(mn.gpts)-num.clust)){
		
		## Clock Expected Computation time
		if(k==1){
			t0 <- proc.time()[1:3]
		}
		
		## Calculate Dissimilarity Between Each Cluster
		if(method!="ward"){
			possible.mergers <- cbind(clust[all.changes[,1]],clust[all.changes[,2]],all.changes[,4])
			diff.clust <- possible.mergers[,1]!=possible.mergers[,2]
			possible.mergers <- possible.mergers[diff.clust,]
			if(method=="complete"){
				possible.mergers <- aggregate(possible.mergers[,3],by=list(clstr1=possible.mergers[,1],clstr2=possible.mergers[,2]),FUN=max)
			} else if(method=="single"){
				possible.mergers <- aggregate(possible.mergers[,3],by=list(clstr1=possible.mergers[,1],clstr2=possible.mergers[,2]),FUN=min)
			} else if(method=="average"){
				possible.mergers <- aggregate(possible.mergers[,3],by=list(clstr1=possible.mergers[,1],clstr2=possible.mergers[,2]),FUN=mean)
			}
			possible.mergers <- possible.mergers[order(possible.mergers[,1]),]
			possible.mergers <- as.matrix(possible.mergers[possible.mergers[,1]<possible.mergers[,2],])
		
		} else {
		
			## Calculate Ward Dissimilarity Between Each Cluster
			cmns <- aggregate(mn.gpts,by=list(cluster=clust),FUN=mean)
			possible.mergers <- aggregate(all.changes[,3],by=list(c1=clust[all.changes[,1]],c2=clust[all.changes[,2]]),FUN=mean)
			possible.mergers <- as.matrix(possible.mergers[possible.mergers[,1]<possible.mergers[,2],])
			clust.sizes <- table(clust)
			ward <- ((cmns[possible.mergers[,2],2]-cmns[possible.mergers[,1],2])^2)*(clust.sizes[possible.mergers[,1]]*clust.sizes[possible.mergers[,2]])/(clust.sizes[possible.mergers[,1]]+clust.sizes[possible.mergers[,2]])
			possible.mergers <- as.matrix(cbind(possible.mergers[,1:2],ward/possible.mergers[,3]))
		
		}
		
		## Determine best merger = smallest dissimilarity
		if(!is.null(dim(possible.mergers))){
			merge.pts <- possible.mergers[possible.mergers[,3]==min(possible.mergers[,3]),1:2]
		} else {
			merge.pts <- possible.mergers[possible.mergers[3]==min(possible.mergers[3])]
		}
    if(length(merge.pts)>2){
      merge.pts <- merge.pts[1,]
    }
	
		## Merge and Recalculate number of clusters
		new.lab <- min(c(merge.pts[1],merge.pts[2]))
#     if(max(c(merge.pts[1],merge.pts[2]))==144){
#       break
#     }
		clust[clust==merge.pts[1]] <- new.lab
		clust[clust==merge.pts[2]] <- new.lab
		clust[clust>max(merge.pts)] <- clust[clust>max(merge.pts)]-1
		
		## Print out expected computation time
		if(k==1){
			t1 <- proc.time()[1:3]-t0
			cat('estimated clustering time (min):',t1[3]*(length(mn.gpts)-num.clust)/60,'...\n')
		}
	} # End k loop
	clust <- as.numeric(factor(clust))
	clust <- clust[assigned.gpt]
	
	cat('Clustering complete.\n')
	
	return(clust)
	
} #End grid.agg.gradclust
	
adj.rand.index <- function(c1,c2) {
  n <- length(c1)
  if ( length(c2) != n ) stop("Clustering must be the same length.")
  t1 <- table(c1)
  t2 <- table(c2)
  t12 <- table(c1,c2)
  expected <- sum(choose(t1,2)) * sum(choose(t2,2)) / choose(n,2)
  numerator <- sum(choose(t12,2)) - expected
  denominator <- 0.5*(sum(choose(t1,2))+sum(choose(t2,2))) - expected
  numerator / denominator
}