## Record the known search limitation; this is not a minimum-size pass.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script), "simulation_support.R"))
source_package_checkout()
args <- parse_args()
output <- arg_value(args,"output",file.path(repo_root(),"simulations/corpus/moment-v2"))
ensure_dir(output)
parameters <- list(k=.1, usd=1, psd=.3, dpm=.03, dpsd=.03, alternative="greater")
curve <- do.call(pnmbf01,c(list(n=1:1000),parameters))
first <- which(curve>=.016)[1]
fast <- do.call(nnmbf01,c(list(power=.016),parameters))
exhaustive <- nbf01seq(k1=.1,k0=10,power=.016,usd=1,psd=.3,dpm=.03,dpsd=.03,
    type="moment",alternative="greater",looks=1,search="exhaustive",nrange=c(2,1000))
result <- data.frame(method=c("nnmbf01","first qualifying integer","nbf01seq exhaustive"),
    n=c(fast,first,exhaustive), power=curve[c(fast,first,exhaustive)],
    minimum_correct=c(fast==first,TRUE,exhaustive==first))
write.csv(result,file.path(output,"known-search-limitation.csv"),row.names=FALSE)
print(result)
stopifnot(first==18,exhaustive==first,curve[first]>=.016,curve[first-1]<.016)
