#################################################################################
# MCWBRMDs_equalsize:Minimal Circular Weakly Balance Repeated Measure Designs 
# for period of equal size(P)

# Algorithm from paper:

#  H.M.Kashif Rasheed, Muhammad Riaz, Muhammad Sajid, Mahmood Ul Hassan, Zahra Noreen
#  and Rashid Ahmed (2021). Algorithm to generate minimal CWBRMDs  
#
# Coded by kashif et al., 2021-2022 
# Version 2.0  (2021-07-25)
#################################################################################




################################################################
# Division of adjusted A in i groups to get the set(s) of shifts
################################################################
grouping1<-function(A,p,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=i){
    
    for(y in 1:5000){
      comp<-sample(1:length(A1),p)
      com<-A1[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs<-rbind(bs,com)
        A1<-A1[-comp]
        z<-z+1
        f=f+1
      }
      if(z==i) break
    }
    if(z<i) {bs<-c();z=0;f=1;A1=A}  
    
  }
  
  
  bs1<-t(apply(bs,1,sort))
  bs1<-cbind(bs1,rowSums(bs),rowSums(bs)/v)
  rownames(bs1)<-paste("G",1:i, sep="")
  colnames(bs1)<-c(paste(1:p, sep=""),"sum" ,"sum/v")
  
  bs2<-t(apply(bs,1,sort))
  bs2<-delmin(bs2)
  list(B1=list(bs2),B2=list(bs1),B3=A1)
}


#######################################################################
# Obtaing set(s) of shifts by deleting smallest value of each group
#######################################################################

delmin<-function(z){
  fs<-c()
  n<-nrow(z)
  c<-ncol(z)-1
  for(i in 1:n){
    z1<-z[i,]
    z2<-z1[z1!=min(z1)]
    fs<-rbind(fs,z2)
  }
  rownames(fs)<-paste("S",1:n, sep="")
  colnames(fs)<-rep("",c)
  return(fs)
}


#################################################################################
# Selection of adjusted A and the set(s) of shifs to obtain Minimal Circular Weakly  
#  Balance RMDs for period of equal size. 
#################################################################################

# D=1: Minimal CGSBRMDs in which v/2 of the ordered pairs do not appear as preceded treatments
# D=2: Minimal CGSBRMDs in which v/2 of the ordered pairs appear twice as preceded treatments
#   P: Period size
#   i: Number of set of shifts for P


CGSBRMD_equalsize<-function(v,p,i,D=1){
  
  if(p<=2) stop("P= Period size: Period size must be greater than 2")
  if(i<=0) stop("i= Must be a positive integer")
  if(v%%2!=0) stop("v must be even")
  
  setClass( "stat_test", representation("list"))
  
  setMethod("show", "stat_test", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following are required sets of shifts to obtain the 
        minimal CWBRMDs for", "v=" ,object[[3]][1], "and","P=",object[[3]][2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    print(object$S[[1]])
  })
  
  if(D==1){  
    
    v=p*i+2
    A<-c(1:((v-2)/2),((v+2)/2),((v+4)/2):(v-1))
    A1<-grouping1(A,p,v,i)
    A2<-c(v,p);names(A2)<-c("V","P")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }
  
  
  if(D==2){
    v=p*i
    A<-c(1:(v-1),(v/2))
    A1<-grouping1(A,p,v,i)
    A2<-c(v,p);names(A2)<-c("V","p")
    x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
  }

  new("stat_test", x)
  
}

##################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from MCWBRMD_equalsize
# The output is called using the design_MCWBRMD to generate design

design_CGSBRMD<-function(H){
  
  setClass( "MCwBRMD_design", representation("list"))
  setMethod("show", "MCwBRMD_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal MCwBRMD for", "v=" ,object$R[1], "and","p=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  })  
  
  v<-H$R[1]
  p<-H$R[2]
  ss<-H$S  
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  
  
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
    
  }
  
  x<-list(Design=G,R=H$R)
  new("CGSBRMD_design", x)
}



##################################################################################
# Examples: Using CSBGND_equalsize function to obtain the set(s) of shifts
# for construction of Minimal Circular weakly Balance Repeated Measure 
# Design for equal period sizes (p)
##################################################################################


# example#
(H<-CGSBRMD_equalsize(v=12,p=5,i=2,D=1))
(D<-design_CGSBRMD(H))

# example #2
(H<-CGSBRMD_equalsize(v=10,p=5,i=2,D=2))
design_CGSBRMD(H)

# example #3
(H<-CGSBRMD_equalsize(v=12,p=3,i=4,D=2))
(H<-CGSBRMD_equalsize(v=6,p=3,i=2,D=2))
H$G
(H<-CGSBRMD_equalsize(p=4,i=2,D=1))
(H<-CGSBRMD_equalsize(p=3,i=2,D=1))
H$G


# example #4
(H<-CGSBRMD_equalsize(p=8,i=7,D=2))
design_CGSBRMD(H)















