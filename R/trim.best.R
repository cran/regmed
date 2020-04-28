trim.best <-
function(obj, lambda=0, mediator.epsilon=0.0001){
    
    if(class(obj) != "regmed.grid") {
      stop("input not regmed.grid class")
    }
	fit.best<-getFit.regmed.grid(obj)

	keep.med<-(abs(fit.best$alpha*fit.best$beta) >= mediator.epsilon)
	keep.med<-paste(paste0("\"",rownames(keep.med)[keep.med],"\""),collapse=",")

	out<-eval(parse(text=paste0("update(fit.best,lambda=",lambda,",mediator=",as.character(fit.best$call)[names(fit.best$call)=="mediator"],"[,c(",keep.med,")])")))

	return(out)
}
