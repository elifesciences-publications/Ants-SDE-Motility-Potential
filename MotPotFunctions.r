
################################################################
##
##   Subroutines
##
################################################################



get.motpot.data <- function(ani.list,boundary.poly,res=1,xcolname="x",ycolname="y",tcolname="t"){
    ######################################
    ## Get raster and adjacency objects
    ######################################
    R=raster(xmn=min(boundary.poly[,1]),xmx=max(boundary.poly[,1]),ymn=min(boundary.poly[,2]),ymx=max(boundary.poly[,2]),res=res)
    values(R) <- NA
    xy=xyFromCell(R,1:ncell(R))
    idx.in.rast=inpip(xy,boundary.poly)
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





motpot.estim <- function(loglam,Q=Q,R=R,ts.df=ts.df,A=A,holdout.idx=integer(),idx.in.rast,cl=NA){
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
        if(is.na(cl)[1]){
            ef=bam(log((as.numeric(eps.hat))^2)~s(x,y),data=ts.df[train.idx,])
        }else{
            ef=bam(log((as.numeric(eps.hat))^2)~s(x,y),data=ts.df[train.idx,],cluster=cl)
        }
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



