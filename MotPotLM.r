################################################################
##
##   Subroutines
##
################################################################

library(Matrix)
library(raster)
library(splancs)


get.motpot.data <- function(ani.list,nest.poly,res=1,xcolname="x",ycolname="y",tcolname="t"){
    ######################################
    ## Get raster and adjacency objects
    ######################################
    R=raster(xmn=min(nest.poly[,1]),xmx=max(nest.poly[,1]),ymn=min(nest.poly[,2]),ymx=max(nest.poly[,2]),res=res)
    values(R) <- NA
    xy=xyFromCell(R,1:ncell(R))
    idx.in.rast=inpip(xy,nest.poly)
    values(R)[idx.in.rast] <- 1
    vals=values(R)
    ##    
    adj=adjacent(R,idx.in.rast,target=idx.in.rast)
    ncell=ncell(R)
    Q=Matrix(0,ncell,ncell,sparse=TRUE)
    Q[adj] <- 1
    ##Q=Q+t(Q)
    one=Matrix(1,nrow=ncell,ncol=1)
    diag(Q)=-Q%*%one
    Q=-Q[idx.in.rast,idx.in.rast]
    n.alpha=nrow(Q)
    ######################################
    ## Make data for motility / potential analysis
    ######################################
    ts.df <- data.frame()
    ants.df <- data.frame()
    cat("Individual Ants ")
    for(i in 1:length(ani.list)){
        cat(i," ")
        x=ani.list[[i]][,xcolname]
        y=ani.list[[i]][,ycolname]
        t=ani.list[[i]][,tcolname]
        cam=ani.list[[i]]$cam
        T=length(x)
        h=1
        vx=x[1:(T-2)]-x[2:(T-1)] ## x(t)-x(t+h)
        wx=(x[3:T]-2*x[2:(T-1)]+x[1:(T-2)])/h ## [ x(t+2h)-2x(t+h)+x(t) ]/h 
        vy=y[1:(T-2)]-y[2:(T-1)] ## y(t)-y(t+h)
        wy=(y[3:T]-2*y[2:(T-1)]+y[1:(T-2)])/h ## [ y(t+2h)-2y(t+h)+y(t) ]/h
        x=x[-(T:(T-1))]
        y=y[-(T:(T-1))]
        t=t[-(T:(T-1))]
        ## cam=cam[-(T:(T-1))]
        ## making "ants.df"
        ants.df <- rbind(ants.df,data.frame(id=names(ani.list)[i],t=t,x=x,y=y,vx=vx,vy=vy,wx=wx,wy=wy))
        ## making "ts.df"
        na.idx=which(is.na(x+y+t+vx+vy+wx+wy))
        nomove.idx=which(vx==0 & vy==0)
        x=x[-unique(c(na.idx,nomove.idx))]
        y=y[-unique(c(na.idx,nomove.idx))]
        t=t[-unique(c(na.idx,nomove.idx))]
        wx=wx[-unique(c(na.idx,nomove.idx))]
        vx=vx[-unique(c(na.idx,nomove.idx))]
        wy=wy[-unique(c(na.idx,nomove.idx))]
        vy=vy[-unique(c(na.idx,nomove.idx))]
        ## cam=cam[-unique(c(na.idx,nomove.idx))]
        ##
        ## Get difference matrix for x-directional derivative of H
        ##
        rx=res(R)[1]
        ry=res(R)[2]
        cell.locs=cellFromXY(R,cbind(x,y))
        cell.locs.up=cellFromXY(R,cbind(x,y+ry))
        cell.locs.down=cellFromXY(R,cbind(x,y-ry))
        cell.locs.right=cellFromXY(R,cbind(x+rx,y))
        cell.locs.left=cellFromXY(R,cbind(x-rx,y))
        idx.x.cent.diff=which(vals[cell.locs.right]==1 & vals[cell.locs.left]==1)
        ## Note: A has dims= (time points , total num cells in R)
        Ax=Matrix(0,length(x),ncell,sparse=TRUE)
        Ax[cbind(idx.x.cent.diff,cell.locs.right[idx.x.cent.diff])] <- 1/2/rx
        Ax[cbind(idx.x.cent.diff,cell.locs.left[idx.x.cent.diff])] <- -1/2/rx
        ##
        idx.y.cent.diff=which(vals[cell.locs.up]==1 & vals[cell.locs.down]==1)
        Ay=Matrix(0,length(y),ncell,sparse=TRUE)
        Ay[cbind(idx.y.cent.diff,cell.locs.up[idx.y.cent.diff])] = 1/2/rx
        Ay[cbind(idx.y.cent.diff,cell.locs.down[idx.y.cent.diff])] = -1/2/rx
        ## remove columns of A that are out of the nest polygon
        idx.x=idx.x.cent.diff
        idx.y=idx.y.cent.diff
        Ax=Ax[idx.x,idx.in.rast]
        Ay=Ay[idx.y,idx.in.rast]
        ## compile data frame
        if(length(idx.x>0)){
            ts.df <- rbind(ts.df,data.frame(id=names(ani.list)[i],t=t[c(idx.x,idx.y)],x=x[c(idx.x,idx.y)],y=y[c(idx.x,idx.y)],w=c(wx[idx.x],wy[idx.y]),v=c(vx[idx.x],vy[idx.y]),y.idx=c(rep(0,length(idx.x)),rep(1,length(idx.y)) )))                                           ##,cam=cam[c(idx.x,idx.y)])
            ##
            if(i==1){
                A=rBind(Ax,Ay)
            }else{
                A=rBind(A,Ax,Ay)
            }
        }
    }
    ##
    ## output
    ##
    list(ants.df=ants.df,ts.df=ts.df,A=A,Q=Q,R=R,idx.in.rast=idx.in.rast)
}





estim <- function(loglam,Q=Q,R=R,ts.df=ts.df,A=A,holdout.idx=integer(),idx.in.rast){
    lambda=exp(loglam)
    QQ=Q
    QQplus=cBind(0,QQ)
    QQplus=rBind(0,QQplus)
    out=list()
    mspe=rep(NA,length(lambda))
    ho=0
    if(length(holdout.idx)>0){
        ho=1
        train.idx=(1:nrow(ts.df))[-holdout.idx]
    }else{
        train.idx=1:nrow(ts.df)
    }
    for(i in length(lambda):1){
        cat(i," ")
        ## get beta.hat and alpha.hat
        X=cBind(ts.df$v[train.idx],A[train.idx,])
        ab.hat=solve(t(X)%*%X+lambda[i]*QQplus,t(X)%*%(ts.df$w[train.idx]))
        ## get resids
        eps.hat=ts.df$w[train.idx]-X%*%ab.hat
        ef=bam(log((as.numeric(eps.hat))^2)~s(x,y),data=ts.df[train.idx,],cluster=cl)
        m=exp(fitted(ef))
        ## re-fit with weighted least squares
        X=cBind(ts.df$v[train.idx]/m,A[train.idx,]/m)
        ab.hat=solve(t(X)%*%X+lambda[i]*QQplus,t(X)%*%(ts.df$w[train.idx]/m))
        eps.hat=ts.df$w[train.idx]-m*X%*%ab.hat
        H.hat=R
        values(H.hat)[idx.in.rast] <- -ab.hat[-1]
        mrastvals=predict(ef,newdata=data.frame(x=xyFromCell(R,idx.in.rast)[,1],y=xyFromCell(R,idx.in.rast)[,2]))
        M=R
        values(M)[idx.in.rast] <- mrastvals
        ## predict
        if(ho==1){
            Xpred=cBind(ts.df$v[-train.idx],A[-train.idx,])
            wpred=Xpred%*%ab.hat
            mspe[i]=mean((wpred-ts.df$w[-train.idx])^2)
        }else{
            mspe[i]=mean(eps.hat^2)
        }    
        out[[i]]=list(ab.hat=ab.hat,H.hat=H.hat,M.hat=M)
    }
    list(out=out,mspe=mspe)
}




############################
##
## Read in Nest Polygons
##
############################

load("/home/ephraim/Dropbox/Ants/Tracking Analysis Hanks/Tracking Plots and Data/Data Products/nest.poly.ld.Rdata")
load("/home/ephraim/Dropbox/Ants/Tracking Analysis Hanks/Tracking Plots and Data/Data Products/nest.poly.hd.Rdata")

##
## Estimation preliminaries
##

library(mgcv)
library(parallel)
cl=makeCluster(4)


############################
##
## Colony 1 low
##
############################


load("/home/ephraim/Dropbox/Ants/Tracking Analysis Hanks/Tracking Plots and Data/Data Products/c1l.list.20170608.Rdata")

## data creation
c1l.mp=get.motpot.data(c1l.list,nest.poly.ld)
names(c1l.mp)
save(c1l.mp,file="c1l.mp.20170608.Rdata")

load("c1l.mp.20170608.Rdata")
ts.df=c1l.mp$ts.df
str(ts.df)



## estimation

holdout.idx.c1l=sample(1:nrow(c1l.mp$ts.df),size=50000)
loglamseq=seq(-8,0,by=1)
eee=estim(loglamseq,c1l.mp$Q,c1l.mp$R,c1l.mp$ts.df,c1l.mp$A,holdout.idx=holdout.idx.c1l,idx.in.rast=c1l.mp$idx.in.rast)
plot(eee$mspe)

idx.min=which.min(eee$mspe)
c1l.mp$lam.best=exp(loglamseq[idx.min])
c1l.mp$H.hat=eee$out[[idx.min]]$H.hat
crs(c1l.mp$H.hat) <- crs("+proj=utm")
c1l.mp$M=eee$out[[idx.min]]$M
c1l.mp$ab.hat=eee$out[[idx.min]]$ab.hat

save(c1l.mp,file="c1l.mp.20170608.Rdata")
load("c1l.mp.20170608.Rdata")
## plots
par(mfrow=c(1,2))
plot(c1l.mp$H.hat,main="Col 1 LD: Potential Surface (Force points downhill)")
plot(c1l.mp$M,main="Col 1 LD: Motility Surface (faster movement in higher motility)")

library(rasterVis)
vectorplot(c1l.mp$H.hat)
plot3D(c1l.mp$H.hat)


############################
##
## Colony 1 high
##
############################

load("/home/ephraim/Dropbox/Ants/Tracking Analysis Hanks/Tracking Plots and Data/Data Products/c1h.list.20170608.Rdata")

## data creation
c1h.mp=get.motpot.data(c1h.list,nest.poly.hd)
names(c1h.mp)
save(c1h.mp,file="c1h.mp.20170608.Rdata")

## estimation

holdout.idx.c1h=sample(1:nrow(c1h.mp$ts.df),size=50000)
loglamseq=seq(-8,0,by=1)
eee=estim(loglamseq,c1h.mp$Q,c1h.mp$R,c1h.mp$ts.df,c1h.mp$A,holdout.idx=holdout.idx.c1h,idx.in.rast=c1h.mp$idx.in.rast)
plot(eee$mspe)

idx.min=which.min(eee$mspe)
c1h.mp$lam.best=exp(loglamseq[idx.min])
c1h.mp$H.hat=eee$out[[idx.min]]$H.hat
crs(c1h.mp$H.hat) <- crs("+proj=utm")
c1h.mp$M=eee$out[[idx.min]]$M
c1h.mp$ab.hat=eee$out[[idx.min]]$ab.hat

save(c1h.mp,file="c1h.mp.20170608.Rdata")

## plots
par(mfrow=c(1,2))
plot(c1h.mp$H.hat,main="Col 1 LD: Potential Surface (Force points downhill)")
plot(c1h.mp$M,main="Col 1 LD: Motility Surface (faster movement in higher motility)")

library(rasterVis)
vectorplot(c1h.mp$H.hat)
plot3D(c1h.mp$H.hat)





############################
##
## Colony 2 low
##
############################

load("/home/ephraim/Dropbox/Ants/Tracking Analysis Hanks/Tracking Plots and Data/Data Products/c2l.list.20170608.Rdata")

## data creation
c2l.mp=get.motpot.data(c2l.list,nest.poly.ld)
names(c2l.mp)
save(c2l.mp,file="c2l.mp.20170608.Rdata")

load("c2l.mp.20170608.Rdata")
## estimation

holdout.idx.c2l=sample(1:nrow(c2l.mp$ts.df),size=50000)
loglamseq=seq(-4,5,by=1)
eee=estim(loglamseq,c2l.mp$Q,c2l.mp$R,c2l.mp$ts.df,c2l.mp$A,holdout.idx=holdout.idx.c2l,idx.in.rast=c2l.mp$idx.in.rast)
plot(eee$mspe)

idx.min=which.min(eee$mspe)
c2l.mp$lam.best=exp(loglamseq[idx.min])
c2l.mp$H.hat=eee$out[[idx.min]]$H.hat
crs(c2l.mp$H.hat) <- crs("+proj=utm")
c2l.mp$M=eee$out[[idx.min]]$M
c2l.mp$ab.hat=eee$out[[idx.min]]$ab.hat

save(c2l.mp,file="c2l.mp.20170608.Rdata")

## plots
par(mfrow=c(1,2))
plot(c2l.mp$H.hat,main="Col 2 LD: Potential Surface (Force points downhill)")
plot(c2l.mp$M,main="Col 2 LD: Motility Surface (faster movement in higher motility)")

library(rasterVis)
vectorplot(c2l.mp$H.hat)
plot3D(c2l.mp$H.hat)



############################
##
## Colony 2 high
##
############################

load("/home/ephraim/Dropbox/Ants/Tracking Analysis Hanks/Tracking Plots and Data/Data Products/c2h.list.20170608.Rdata")

## data creation
c2h.mp=get.motpot.data(c2h.list,nest.poly.hd)
names(c2h.mp)
save(c2h.mp,file="c2h.mp.20170608.Rdata")

## estimation

holdout.idx.c2h=sample(1:nrow(c2h.mp$ts.df),size=50000)
loglamseq=seq(-8,0,by=1)
eee=estim(loglamseq,c2h.mp$Q,c2h.mp$R,c2h.mp$ts.df,c2h.mp$A,holdout.idx=holdout.idx.c2h,idx.in.rast=c2h.mp$idx.in.rast)
plot(eee$mspe)

idx.min=which.min(eee$mspe)
c2h.mp$lam.best=exp(loglamseq[idx.min])
c2h.mp$H.hat=eee$out[[idx.min]]$H.hat
crs(c2h.mp$H.hat) <- crs("+proj=utm")
c2h.mp$M=eee$out[[idx.min]]$M
c2h.mp$ab.hat=eee$out[[idx.min]]$ab.hat

save(c2h.mp,file="c2h.mp.20170608.Rdata")

## plots
par(mfrow=c(1,2))
plot(c2h.mp$H.hat,main="Col 2 LD: Potential Surface (Force points downhill)")
plot(c2h.mp$M,main="Col 2 LD: Motility Surface (faster movement in higher motility)")

library(rasterVis)
vectorplot(c2h.mp$H.hat)
plot3D(c2h.mp$H.hat)





############################
##
## Colony 3 low
##
############################

load("/home/ephraim/Dropbox/Ants/Tracking Analysis Hanks/Tracking Plots and Data/Data Products/c3l.list.20170608.Rdata")

## data creation
c3l.mp=get.motpot.data(c3l.list,nest.poly.ld)
names(c3l.mp)
save(c3l.mp,file="c3l.mp.20170608.Rdata")

## estimation

holdout.idx.c3l=sample(1:nrow(c3l.mp$ts.df),size=50000)
loglamseq=seq(-8,0,by=1)
eee=estim(loglamseq,c3l.mp$Q,c3l.mp$R,c3l.mp$ts.df,c3l.mp$A,holdout.idx=holdout.idx.c3l,idx.in.rast=c3l.mp$idx.in.rast)
plot(eee$mspe)

idx.min=which.min(eee$mspe)
c3l.mp$lam.best=exp(loglamseq[idx.min])
c3l.mp$H.hat=eee$out[[idx.min]]$H.hat
crs(c3l.mp$H.hat) <- crs("+proj=utm")
c3l.mp$M=eee$out[[idx.min]]$M
c3l.mp$ab.hat=eee$out[[idx.min]]$ab.hat

save(c3l.mp,file="c3l.mp.20170608.Rdata")

## plots
par(mfrow=c(1,2))
plot(c3l.mp$H.hat,main="Col 3 LD: Potential Surface (Force points downhill)")
plot(c3l.mp$M,main="Col 3 LD: Motility Surface (faster movement in higher motility)")

library(rasterVis)
vectorplot(c3l.mp$H.hat)
plot3D(c3l.mp$H.hat)




############################
##
## Colony 3 high
##
############################

load("/home/ephraim/Dropbox/Ants/Tracking Analysis Hanks/Tracking Plots and Data/Data Products/c3h.list.20170608.Rdata")

## data creation
c3h.mp=get.motpot.data(c3h.list,nest.poly.hd)
names(c3h.mp)
save(c3h.mp,file="c3h.mp.20170608.Rdata")

## estimation

holdout.idx.c3h=sample(1:nrow(c3h.mp$ts.df),size=50000)
loglamseq=seq(-8,0,by=1)
eee=estim(loglamseq,c3h.mp$Q,c3h.mp$R,c3h.mp$ts.df,c3h.mp$A,holdout.idx=holdout.idx.c3h,idx.in.rast=c3h.mp$idx.in.rast)
plot(eee$mspe)

idx.min=which.min(eee$mspe)
c3h.mp$lam.best=exp(loglamseq[idx.min])
c3h.mp$H.hat=eee$out[[idx.min]]$H.hat
crs(c3h.mp$H.hat) <- crs("+proj=utm")
c3h.mp$M=eee$out[[idx.min]]$M
c3h.mp$ab.hat=eee$out[[idx.min]]$ab.hat

save(c3h.mp,file="c3h.mp.20170608.Rdata")

## plots
par(mfrow=c(1,2))
plot(c3h.mp$H.hat,main="Col 3 LD: Potential Surface (Force points downhill)")
plot(c3h.mp$M,main="Col 3 LD: Motility Surface (faster movement in higher motility)")

library(rasterVis)
vectorplot(c3h.mp$H.hat)
plot3D(c3h.mp$H.hat)





################################3
##
## plots for all colonies
##
################################


load("c1h.mp.20170608.Rdata")
load("c1l.mp.20170608.Rdata")
load("c2h.mp.20170608.Rdata")
load("c2l.mp.20170608.Rdata")
load("c3h.mp.20170608.Rdata")
load("c3l.mp.20170608.Rdata")

library(fields)
library(rasterVis)
library(RColorBrewer)


?colorRampPalette

rb.colors=colorRampPalette(c("red","blue"))
yp.colors=colorRampPalette(c("yellow","purple"))
##
## Motility Plots
##

pdf("C1mot.pdf",width=20,height=5)
layout(matrix(c(1,2,2,2,2),1))
image( c1h.mp$M, col=( yp.colors( 99 ) ), breaks=seq(min(minValue( c1l.mp$M )),max(maxValue(c1l.mp$M)),length.out=100), bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C1 High Density Motility") 
##
image.plot( c1l.mp$M, col=( yp.colors( 99 ) ), breaks=seq(min(minValue( c1l.mp$M )),max(maxValue(c2l.mp$M)),length.out=100) ,bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C1 Low Density Motility") 
dev.off()

pdf("C2mot.pdf",width=20,height=5)
layout(matrix(c(1,2,2,2,2),1))
image( c2h.mp$M, col=( yp.colors( 99 ) ), breaks=seq(min(minValue( c2h.mp$M )),max(maxValue(c2l.mp$M)),length.out=100), bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C2 High Density Motility") 
##
image.plot( c2l.mp$M, col=( yp.colors( 99 ) ), breaks=seq(min(minValue( c2l.mp$M )),max(maxValue(c2l.mp$M)),length.out=100) ,bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C2 Low Density Motility") 
dev.off()


pdf("C3mot.pdf",width=20,height=5)
layout(matrix(c(1,2,2,2,2),1))
image( c3h.mp$M, col=( yp.colors( 99 ) ), breaks=seq(min(minValue( c2l.mp$M )),max(maxValue(c2l.mp$M)),length.out=100), bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C3 High Density Motility") 
##
image.plot( c3l.mp$M, col=( yp.colors( 99 ) ), breaks=seq(min(minValue( c2l.mp$M )),max(maxValue(c2l.mp$M)),length.out=100) ,bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C3 Low Density Motility") 
dev.off()


##
## Potential surface plots
##

## centering all potential surfaces
c1h.mp$H.hat=c1h.mp$H.hat-min(values(c1h.mp$H.hat),na.rm=TRUE)
c1l.mp$H.hat=c1l.mp$H.hat-min(values(c1l.mp$H.hat),na.rm=TRUE)
c2h.mp$H.hat=c2h.mp$H.hat-min(values(c2h.mp$H.hat),na.rm=TRUE)
c2l.mp$H.hat=c2l.mp$H.hat-min(values(c2l.mp$H.hat),na.rm=TRUE)
c3h.mp$H.hat=c3h.mp$H.hat-min(values(c3h.mp$H.hat),na.rm=TRUE)
c3l.mp$H.hat=c3l.mp$H.hat-min(values(c3l.mp$H.hat),na.rm=TRUE)

## finding max
max.val=max(values(c1h.mp$H.hat),na.rm=TRUE)
max.val=max(c(max.val,values(c1l.mp$H.hat)),na.rm=TRUE)
max.val=max(c(max.val,values(c2l.mp$H.hat)),na.rm=TRUE)
max.val=max(c(max.val,values(c2h.mp$H.hat)),na.rm=TRUE)
max.val=max(c(max.val,values(c3l.mp$H.hat)),na.rm=TRUE)
max.val=max(c(max.val,values(c3h.mp$H.hat)),na.rm=TRUE)
max.val

## centering potentials for high density
c1h.mp$H.hat=c1h.mp$H.hat-mean(values(c1h.mp$H.hat),na.rm=TRUE)+65
c2h.mp$H.hat=c2h.mp$H.hat-mean(values(c2h.mp$H.hat),na.rm=TRUE)+65
c3h.mp$H.hat=c3h.mp$H.hat-mean(values(c3h.mp$H.hat),na.rm=TRUE)+65


library(ctmcmove)

quiver <- function(rast,spacing=1,scaling=1,...){
    rast=aggregate(rast,spacing)
    R=rast.grad(rast)
    R2=rast
    ##R2=aggregate(rast,spacing)
    xy=xyFromCell(R2,cell=1:ncell(R2))
    cells=cellFromXY(rast,xy)
    idx.in=which(values(rast)[cells]>-Inf)
    cells=cells[idx.in]
    xy.start=xy[idx.in,]
    xy.end=xy.start
    xy.end[,1]=xy.end[,1]+R$rast.grad.x[cells]*scaling
    xy.end[,2]=xy.end[,2]+R$rast.grad.y[cells]*scaling
    length.arrow=apply(abs(xy.end-xy.start),1,sum)
    idx.nonzero=which(length.arrow>0)
    arrows(xy.start[idx.nonzero,1],xy.start[idx.nonzero,2],xy.end[idx.nonzero,1],xy.end[idx.nonzero,2],...)
}




## figure supplement 3

image( c1h.mp$H.hat, col=( yp.colors( 99 ) ), breaks=seq(0,max.val,length.out=100), bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C1 High Density Potential") 
quiver(c1h.mp$H.hat,spacing=2,scaling=-3,length=.04,lwd=.6)
##
image.plot( c1l.mp$H.hat, col=( yp.colors( 99 ) ), breaks=seq(0,max.val,length.out=100) ,bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C1 Low Density Potential") 
quiver(c1l.mp$H.hat,spacing=2,scaling=-3,length=.04,lwd=.6)



pdf("C1pot.pdf",width=20,height=5)
image( c1h.mp$H.hat, col=( yp.colors( 99 ) ), breaks=seq(0,max.val,length.out=100), bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C1 High Density Potential") 
quiver(c1h.mp$H.hat,spacing=2,scaling=-3,length=.04,lwd=.6)
##
image.plot( c1l.mp$H.hat, col=( yp.colors( 99 ) ), breaks=seq(0,max.val,length.out=100) ,bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C1 Low Density Potential") 
quiver(c1l.mp$H.hat,spacing=2,scaling=-3,length=.04,lwd=.6)
dev.off()



pdf("C2pot.pdf",width=20,height=5)
image( c2h.mp$H.hat, col=( yp.colors( 99 ) ), breaks=seq(0,max.val,length.out=100), bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C2 High Density Potential") 
quiver(c2h.mp$H.hat,spacing=2,scaling=-3,length=.04,lwd=.6)
##
image.plot( c2l.mp$H.hat, col=( yp.colors( 99 ) ), breaks=seq(0,max.val,length.out=100) ,bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C2 Low Density Potential") 
quiver(c2l.mp$H.hat,spacing=2,scaling=-3,length=.04,lwd=.6)
dev.off()




pdf("C3pot.pdf",width=20,height=5)
image( c3h.mp$H.hat, col=( yp.colors( 99 ) ), breaks=seq(0,max.val,length.out=100), bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C3 High Density Potential") 
quiver(c3h.mp$H.hat,spacing=2,scaling=-3,length=.04,lwd=.6)
##
image.plot( c3l.mp$H.hat, col=( yp.colors( 99 ) ), breaks=seq(0,max.val,length.out=100) ,bty="n",xlab="",ylab="",axes=FALSE,asp=1 ,main="C3 Low Density Potential") 
quiver(c3l.mp$H.hat,spacing=2,scaling=-3,length=.04,lwd=.6)
dev.off()









##
## Data Plots
##

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col.vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col.vector, n))

col.vector=c(3,7,8,12,17,26,30,31,32,33,
            40,41,45,47,51,53,56,57,76,78,
            81,84,96,100,101,102,
            367,370,371,372,373,375,376,382,399,419,
            426,428,429,430,434,435,450,454,454,461,
            462,466,467,471,476,482,490,493,498,501,
            502,503,506,507,508,512,513,514,518,519,
            551,552,556,563,567,571,574,575,610,614,615,616,
            620,625,630,640,645,650,655,660,665,670,675,680)

## colony 1


png("c1l.png",width=15000,height=12000)#,width=8.5,height=11)
par(mfrow=c(2,1))
plot(nest.poly.ld,type="l",asp=1,xlab="mm",ylab="mm",lwd=50)
for(i in 1:length(c1l.list)){
    points(c1l.list[[i]][-1,2:3],col=col.vector[i],pch=20,cex=4,type="p")
}
dev.off()

pdf("c1l.pdf")#,width=8.5,height=11)
par(mfrow=c(2,1))
for(i in 1:length(c1l.list)){
    plot(nest.poly.ld,type="l",asp=1,main=names(c1l.list)[i],xlab="mm",ylab="mm")
    points(c1l.list[[i]][-1,2:3],col=col.vector[i],pch=".",type="p",cex=2)
    points(c1l.list[[i]][-1,2:3],col=col.vector[i],type="l",cex=2)
}
dev.off()


png("c1h.png",width=15000,height=12000)#,width=8.5,height=11)
par(mfrow=c(2,1))
plot(nest.poly.hd,type="l",asp=1,xlab="mm",ylab="mm",lwd=50)
for(i in 1:length(c1h.list)){
    points(c1h.list[[i]][-1,2:3],col=col.vector[i],pch=20,cex=4,type="p")
}
dev.off()

pdf("c1h.pdf")#,width=8.5,height=11)
par(mfrow=c(2,2))
for(i in 1:length(c1h.list)){
    plot(nest.poly.hd,type="l",asp=1,main=names(c1h.list)[i],xlab="mm",ylab="mm")
    points(c1h.list[[i]][-1,2:3],col=col.vector[i],pch=".",type="p",cex=2)
    points(c1h.list[[i]][-1,2:3],col=col.vector[i],type="l",cex=2)
}
dev.off()



## colony 2



png("c2l.png",width=15000,height=12000)#,width=8.5,height=11)
par(mfrow=c(2,1))
plot(nest.poly.ld,type="l",asp=1,xlab="mm",ylab="mm",lwd=50)
for(i in 1:length(c2l.list)){
    points(c2l.list[[i]][-1,2:3],col=col.vector[i],pch=20,cex=4,type="p")
}
dev.off()

pdf("c2l.pdf")#,width=8.5,height=11)
par(mfrow=c(2,1))
for(i in 1:length(c2l.list)){
    plot(nest.poly.ld,type="l",asp=1,main=names(c2l.list)[i],xlab="mm",ylab="mm")
    points(c2l.list[[i]][-1,2:3],col=col.vector[i],pch=".",type="p",cex=2)
    points(c2l.list[[i]][-1,2:3],col=col.vector[i],type="l",cex=2)
}
dev.off()


png("c2h.png",width=15000,height=12000)#,width=8.5,height=11)
par(mfrow=c(2,1))
plot(nest.poly.hd,type="l",asp=1,xlab="mm",ylab="mm",lwd=50)
for(i in 1:length(c2h.list)){
    points(c2h.list[[i]][-1,2:3],col=col.vector[i],pch=20,cex=4,type="p")
}
dev.off()

pdf("c2h.pdf")#,width=8.5,height=11)
par(mfrow=c(2,2))
for(i in 1:length(c2h.list)){
    plot(nest.poly.hd,type="l",asp=1,main=names(c2h.list)[i],xlab="mm",ylab="mm")
    points(c2h.list[[i]][-1,2:3],col=col.vector[i],pch=".",type="p",cex=2)
    points(c2h.list[[i]][-1,2:3],col=col.vector[i],type="l",cex=2)
}
dev.off()



## colony 3

png("c3l.png",width=15000,height=12000)#,width=8.5,height=11)
par(mfrow=c(2,1))
plot(nest.poly.ld,type="l",asp=1,xlab="mm",ylab="mm",lwd=50)
for(i in 1:length(c3l.list)){
    points(c3l.list[[i]][-1,2:3],col=col.vector[i],pch=20,cex=4,type="p")
}
dev.off()

pdf("c3l.pdf")#,width=8.5,height=11)
par(mfrow=c(2,1))
for(i in 1:length(c3l.list)){
    plot(nest.poly.ld,type="l",asp=1,main=names(c3l.list)[i],xlab="mm",ylab="mm")
    points(c3l.list[[i]][-1,2:3],col=col.vector[i],pch=".",type="p",cex=2)
    points(c3l.list[[i]][-1,2:3],col=col.vector[i],type="l",cex=2)
}
dev.off()


png("c3h.png",width=15000,height=12000)#,width=8.5,height=11)
par(mfrow=c(2,1))
plot(nest.poly.hd,type="l",asp=1,xlab="mm",ylab="mm",lwd=50)
for(i in 1:length(c3h.list)){
    points(c3h.list[[i]][-1,2:3],col=col.vector[i],pch=20,cex=4,type="p")
}
dev.off()

pdf("c3h.pdf")#,width=8.5,height=11)
par(mfrow=c(2,2))
for(i in 1:length(c3h.list)){
    plot(nest.poly.hd,type="l",asp=1,main=names(c3h.list)[i],xlab="mm",ylab="mm")
    points(c3h.list[[i]][-1,2:3],col=col.vector[i],pch=".",type="p",cex=2)
    points(c3h.list[[i]][-1,2:3],col=col.vector[i],type="l",cex=2)
}
dev.off()

