match.nn <- function(A,dt,ps,M=1,both=1){
    if (both==1){
        n <- length(A)
        matchset <- rbindlist(lapply(1:n,function(i){
            ai <- A[i]
            psi <- ps[i]
            dt.i <- dt[A!=ai]
            dps.i <- abs(ps[A!=ai]-psi)
            dt.i[,dps:=dps.i]
            setkey(dt.i,dps)
            out <- rbindlist(list(cbind(dt[i],dps=0),dt.i[1:M,]),use.names=TRUE)
            out[,index:=i]
            out
        }))
        matchset
    }else{
        n1 <- sum(A==1)
        casedt <- dt[A==1]
        controldt <- dt[A==0]
        caseps <- ps[A==1]
        controlps <- ps[A==0]
        casedt[,ps:=caseps]
        controldt[,ps:=controlps]
        matchset <- rbindlist(lapply(1:n1,function(i){
            casei <- casedt[i]
            psi <- caseps[i]
            dpsi <- abs(controlps-psi)
            matchi <- controldt[order(dpsi)][1:M]
            out <- rbind(casei,matchi)
                        out[,index:=i]
            out
        }))
        matchset
    }
}

match.fun <- function(A,Y,ps,M=5){
    n <- length(A)
    dt <- data.table(Y,A,ps)
    Yhat <- sapply(1:n,function(i){
        ai <- A[i]
        psi <- ps[i]
        dt.i <- dt[A!=ai]
        dps.i <- abs(ps[A!=ai]-psi)
        dt.i[,dps:=dps.i]
        setkey(dt.i,dps)
        # match.set <- dt.i[1:M]
        Yhat.i <- dt.i[1:M,mean(Y)]
        ## if (i%in%c(7,78,348,1322,2214)) {
            ## print(paste0("Case ",i," A[i]=",A[i]," Y[i]=",Y[i]," Yihat=",Yhat.i," ps[i]=",psi))
            ## print(dt.i[1:M])
        ## }
        Yhat.i
    })
    1/n * sum((2*A-1)*(Y-Yhat))
}
sd.match.fun <- function(A,Y,ps,theta,X,M=5,L=2){
    n <- length(A)
    eta <- theta[1]+X%*%theta[-1]
    F.fun <- 1/(1+exp(-eta))
    f.fun <- exp(-eta)/(1+exp(-eta))^2
    I.hat <- 1/n * Reduce("+",lapply(1:n,function(i){
        (f.fun^2/(F.fun*(1-F.fun)))[i]* X[i,]%*%t(X[i,])
    }))
    HL.i.set <- lapply(1:n,function(i){
        ai <- A[i]
        in.ai <- (1:n)[A==ai]
        tmp <- data.frame(id.ai=in.ai,X=X[in.ai,],Y=Y[in.ai],ps=ps[in.ai],d.ps.ai=abs(ps[in.ai]-ps[i]))
        return(tmp[order(tmp$d.ps.ai),][1:L,])
    })
    JM.i.set <- lapply(1:n,function(i){
        ai <- A[i]
        not.ai <- (1:n)[A!=ai]
        tmp <- data.frame(id.not.ai=not.ai,X=X[not.ai,],Y=Y[not.ai],ps=ps[not.ai],d.ps.ai=abs(ps[not.ai]-ps[i]))
        return(tmp[order(tmp$d.ps.ai),][1:M,])
    })
    JL.i.set <- lapply(1:n,function(i){
        ai <- A[i]
        not.ai <- (1:n)[A!=ai]
        tmp <- data.frame(id.not.ai=not.ai,X=X[not.ai,],Y=Y[not.ai],ps=ps[not.ai],d.ps.ai=abs(ps[not.ai]-ps[i]))
        return(tmp[order(tmp$d.ps.ai),][1:L,])
    })
    sigma.bar.hat <- 1/(L-1)*sapply(1:n,function(i){
        sum((HL.i.set[[i]]$Y-mean(HL.i.set[[i]]$Y))^2)
    })
    tau.hat <- 1/n * sum((2*A-1)*(Y-1/M * sapply(1:n,function(i){sum(JM.i.set[[i]]$Y)})))
    JM.full.set <- unlist(lapply(1:n,function(i){
        JM.i.set[[i]]$id.not.ai
    }))
    KM <- sapply(1:n,function(j){sum(JM.full.set==j)})
    sigma.hat <- 1/n * sum(((2*A-1)*(Y-1/M*sapply(1:n,function(i){sum(JM.i.set[[i]]$Y)}))-tau.hat)^2)+
        1/n * sum( ((KM/M)^2+((2*M-1)/M)*(KM/M))*sigma.bar.hat)
    cov.W.w <- 1/(L-1)*do.call(rbind,lapply(1:n,function(i){
        HL.tmp <- HL.i.set[[i]]
        X.i <- HL.tmp[,substr(names(HL.tmp),1,1)=="X"]
        return(sapply(1:ncol(X.i),function(j){
            sum((X.i[,j]-mean(X.i[,j]))*(HL.tmp$Y-mean(HL.tmp$Y)))
        }))
    }))
    cov.W.not.w <- 1/(L-1)*do.call(rbind,lapply(1:n,function(i){
        JL.tmp <- JL.i.set[[i]]
        X.i <- JL.tmp[,substr(names(JL.tmp),1,1)=="X"]
        return(sapply(1:ncol(X.i),function(j){
            sum((X.i[,j]-mean(X.i[,j]))*(JL.tmp$Y-mean(JL.tmp$Y)))
        }))
    }))
    chat <- 1/n* sapply(1:ncol(X),function(j){
        sum((((A==1)*cov.W.w+(A==0)*cov.W.not.w)/ps+((A==0)*cov.W.w+(A==1)*cov.W.not.w)/(1-ps))[,j]*matrix(f.fun,1,n)
            )})
    final.sigma <- sigma.hat-matrix(chat,1,ncol(X))%*%solve(I.hat)%*%matrix(chat,ncol(X),1)
    final.sigma <- 1/sqrt(n)*sqrt(final.sigma[1])
    attr(final.sigma,"unadj") <- 1/sqrt(n)*sqrt(sigma.hat)
    return(final.sigma)
}

## *** Test in simulated data
## library(riskRegression)
## set.seed(17)
## d <- sampleData(100,outcome="binary")
## f <- glm(X1~X6+X2,data=d,family=binomial())
## p <- predictRisk(f,newdata=d)
## d$X1 <- as.numeric(d$X1)
## mean(d$X1)
## x <- model.design(terms(model.frame(~X6+X2,data=d)),data=d)
## match.fun(A=d$X1,Y=d$Y,ps=p,M=5)
## sd.match.fun(A=d$X1,Y=d$Y,ps=p,theta=coef(f),X=x$design[,-1],M=5,L=2)
