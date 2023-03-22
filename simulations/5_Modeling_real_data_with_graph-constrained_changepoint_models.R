#(method 1) ON GITHUB
#devtools::install_github("vrunge/gfpop", force = TRUE)

#(method 2) ON CRAN
install.packages("gfpop")
library(gfpop)

devtools::install_github("vrunge/gfpop.data")
library(gfpop.data)
data(profile614chr2, package = "gfpop.data")
data(neuroSpike, package = "gfpop.data")
data(ECG, package = "gfpop.data")

library(penaltyLearning)
library(data.table)
library(ggplot2)

## ----------------------------------------------- ##
## Relevant changes model for DNA copy number data ##
## ----------------------------------------------- ##

profile614chr2$probes[, change.after := floor(
  c(diff(position)/2+position[-.N], NA)) ]
sngraph <- function(n.segs, type, gap){
  stopifnot(is.integer(n.segs), length(n.segs)==1, n.segs >= 1)
  s <- n.segs-1
  seg.vec <- 1:s
  gfpop::graph(
    gfpop::StartEnd(start=0, end=s),
    gfpop::Edge(seg.vec-1, seg.vec, type, gap=gap), all.null.edges = TRUE)
}
select.dt <- rbind(
  data.table(n.segments = c(3, 7, 13), graph.name = "std", gap = 0),
  data.table(n.segments = 13, graph.name = "abs", gap = 1))
seg.dt.list <- list()
for(model.i in 1:nrow(select.dt)){
  model.info <- select.dt[model.i]
  cat(model.i, nrow(select.dt))
  g <- model.info[, sngraph(as.integer(n.segments), graph.name, gap=gap)]
  fit <- gfpop::gfpop(
    profile614chr2$probes$logratio,
    mygraph = g, type = "mean")
  end.i <- fit$changepoints
  change.i <- end.i[-length(end.i)]
  change.pos <- profile614chr2$probes$change.after[change.i]
  seg.dt.list[[model.i]] <- data.table(
    model.info,
    segStart=c(profile614chr2$probes[1, position], change.pos),
    segEnd=c(change.pos, profile614chr2$probes[.N, position]),
    mean=fit$parameters)
}
seg.dt <- do.call(rbind, seg.dt.list)

some.segs <- seg.dt[select.dt, on=list(n.segments, graph.name), allow.cartesian=TRUE]
some.change <- some.segs[min(segStart) < segStart]
some.change[, change := segStart]

err.dt <- some.segs[, {
  change.dt <- .SD[min(segStart) < segStart]
  change.dt[, change := segStart]
  change.dt[, prob := 1]
  change.dt[, nseg := n.segments]
  model.dt <- data.table(nseg=n.segments, prob=1)
  penaltyLearning::labelError(
    model.dt,
    data.table(prob=1, profile614chr2$labels),
    change.dt,
    change.var="change",
    model.vars="nseg",
    label.vars=c("labelStart", "labelEnd"),
    problem.vars="prob")$label.errors
}, by=list(graph.name, n.segments)]

err.dt[, list(
  fp=sum(fp),
  fn=sum(fn),
  errors=sum(fp+fn)
  ), by=list(graph.name, n.segments)]

win <- function(min,max)data.table(windowStart=min*1e5,windowEnd=max*1e5)
windows <- rbind(
  win(  65,  71),
  win( 148, 171),
  win( 354, 361),
  win(1059,1065))
mb.fac <- 1e6
wfac <- function(x){
  factor(x, c("6.5-7.1", "14.8-17.1", "35.4-36.1", "105.9-106.5"))
}
windows[, window := wfac(windowStart/mb.fac, windowEnd/mb.fac)]
setkey(windows, windowStart, windowEnd)
f <- function(dt, key.vec){
  setkeyv(dt, key.vec)
  dt
}
profile614chr2$probes[, pos0 := position]
over.list <- list(
  changes=f(some.change, c("change", "segStart")),
  segments=f(some.segs, c("segStart", "segEnd")),
  labels=f(profile614chr2$labels, c("labelStart", "labelEnd")),
  errors=f(err.dt, c("labelStart", "labelEnd")),
  probes=f(profile614chr2$probes, c("position", "pos0")))
join.list <- lapply(over.list, foverlaps, windows, nomatch=0L)

show.err <- join.list$errors[, list(
  fp=sum(fp),
  fn=sum(fn),
  errors=sum(fp+fn)
  ), by=list(graph.name, n.segments)]

br <- c(6.5,7.0,seq(15,17,by=0.5),35.5,36,106,106.5)
names(br) <- as.character(br)
gg.out <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  coord_cartesian(expand=FALSE)+
  facet_grid(
    graph.name + n.segments ~ window,
    ##labeller=label_both,
    scales="free",space="free")+
  penaltyLearning::geom_tallrect(aes(
    xmin=labelStart/mb.fac, xmax=labelEnd/mb.fac, fill=annotation),
    color="grey",
    data=join.list$labels)+
  scale_fill_manual(
    "label",
    values=penaltyLearning::change.colors)+
  penaltyLearning::geom_tallrect(aes(
    xmin=labelStart/mb.fac, xmax=labelEnd/mb.fac, linetype=status),
    fill=NA,
    size=1,
    color="black",
    data=join.list$errors)+
  scale_linetype_manual(
    "error type",
    limits=c("correct",
             "false negative",
             "false positive"),
    values=c(correct=0,
             "false negative"=3,
             "false positive"=1))+
  geom_point(aes(
    position/mb.fac, logratio),
    shape=1,
    color="grey50",
    data=join.list$probes)+
  scale_y_continuous(
    "Logratio (Measure of DNA copy number)",
    breaks=seq(-2, 2, by=2)
  )+
  scale_x_continuous("Position on chr2 (mega bases)", breaks=br)+
  scale_color_manual(values=c(
    std="green",
    abs="deepskyblue"))+
  geom_segment(aes(
    ifelse(segStart < windowStart, windowStart, segStart)/mb.fac, mean,
    color=graph.name,
    xend=ifelse(windowEnd < segEnd, windowEnd, segEnd)/mb.fac, yend=mean),
    data=join.list$segments,
    size=1)+
  geom_vline(aes(
    xintercept=change/mb.fac, color=graph.name),
    linetype="dashed",
    size=0.75,
    data=join.list$changes)+
  geom_text(aes(
    14.8, -3, label=cat(
      " %d label error%s for %d segment %s model",
      errors, ifelse(errors==1, "", "s"), n.segments, graph.name)),
    hjust=0,
    vjust=0,
    data=data.table(show.err, window=wfac("14.8-17.1")))

pdf( "figure11.pdf" )
print(gg.out)
dev.off()

## ----------------------------------------------------- ##
## Multi-modal regression for neuro spike train data set ##
## ----------------------------------------------------- ##

fps <- 100
seconds.between.data <- 1/fps
sec.w <- seconds.between.data/3
myGraph <- gfpop::graph(
  gfpop::Edge(0, 0, "down", penalty = 0),
  gfpop::Edge(1, 1, "up",  penalty =  0),
  gfpop::Edge(0, 1, "up",  penalty =   5),
  gfpop::Edge(1, 0, "down",  penalty = 0))

fit <- gfpop::gfpop(
  data = neuroSpike$calcium, mygraph = myGraph, type = "mean")

end.i <- fit$changepoints
start.i <- c(1, end.i[-length(end.i)]+1)
res.dt <- with(neuroSpike, data.table(
  start.seconds=seconds[start.i],
  end.seconds=seconds[end.i],
  state=fit$states,
  mean=fit$parameters))
res.dt[, Multimodal := ifelse(mean<0.1, 0, mean) ]
over.dt <- res.dt[neuroSpike, on=list(
  start.seconds <= seconds, end.seconds >= seconds)]
tall.dt <- melt(
  over.dt,
  measure.vars=c("Multimodal", "AR1"),
  variable.name="model")
tall.dt[, change.after := c(diff(value), NA)]
tall.dt[, seconds.after := c(start.seconds[-1], NA)]
tall.dt[, spike.i := cumsum(change.after < 0)]
tall.dt[, thresh := ifelse(model=="Multimodal", 0, 0)]
m <- function(model){
  factor(
    model,
    c("AR1", "Multimodal"),
    c("Previous model:
AR1 changepoint
Jewell et al 2017", "Proposed model:
changepoint with
graph constraints"))
}
tall.dt[, model.fac := m(model)]
spike.dt <- tall.dt[0 < change.after & thresh < value, list(
  start.seconds=min(start.seconds),
  end.seconds=max(start.seconds)
), by=list(spike.i, model.fac)]
spike.dt[, mid.seconds := (start.seconds+end.seconds)/2]

lab <- function(xmin, xmax){
  data.table(xmin, xmax, label="oneSpike", annotation="1change", problem=1)
}
label.dt <- rbind(
  lab(166.5, 166.75),
  lab(169.7, 169.9))
ann.colors <- c(oneSpike="#ffafaf")
xmin <- 166
xmax <- 171
spike.dt[, problem := 1]
models.dt <- spike.dt[, list(
  models=.N
  ), by=list(model.fac, problem)]
err.list <- penaltyLearning::labelError(
  models.dt, label.dt, spike.dt,
  change.var="mid.seconds",
  problem.var="problem",
  label.vars=c("xmin", "xmax"),
  model.vars="model.fac")
type.colors <- c(
  data="grey50",
  model="blue")
show.spikes <- spike.dt[xmin < mid.seconds & mid.seconds < xmax]
show.data <- neuroSpike[xmin < seconds & seconds < xmax]
spike.y <- -1.5
gg.out <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(model.fac ~ .)+
  penaltyLearning::geom_tallrect(aes(
    xmin=xmin, xmax=xmax, fill=label),
    color=NA,
    data=label.dt)+
  penaltyLearning::geom_tallrect(aes(
    xmin=xmin, xmax=xmax, linetype=status),
    fill=NA,
    color="black",
    data=err.list$label.errors)+
  scale_linetype_manual(
    "error type",
    limits=c("correct",
             "false negative",
             "false positive"),
    values=c(correct=0,
             "false negative"=3,
             "false positive"=1))+
  geom_point(aes(
    seconds, calcium, color=type),
    shape=1,
    data=data.table(type="data", show.data))+
  geom_line(aes(
    start.seconds, value, color=type),
    data=data.table(
      type="model",
      tall.dt[xmin < start.seconds & start.seconds < xmax]),
    size=0.5)+
  geom_point(aes(
    mid.seconds, spike.y, color=type),
    shape=1,
    data=data.table(type="model", show.spikes))+
  scale_y_continuous(
    "Fluorescence intensity
(Measure of neural activity)",
    breaks=seq(0, 10, by=2),
    limits=c(-2, 10)
  )+
  scale_x_continuous("Time (seconds)")+
  guides(color="none")+
  scale_color_manual(values=type.colors)+
  geom_text(aes(
    x,y,label=label, color=type, hjust=hjust),
    size=3,
    data=data.table(model.fac=m("AR1"), rbind(
      data.table(
        hjust=0, x=167, y=2, label="Noisy activity data", type="data"),
      data.table(
        hjust=1, x=169.5, y=5, label="Mean model", type="model"),
      data.table(
        hjust=1, x=169.4, y=spike.y, label="Predicted spikes", type="model"))))

pdf( "figure12.pdf" )
print(gg.out)
dev.off()


## -------------------------------------------------------------------------- ##
## Nine-state model for QRS complex detection in electrocardiogram (ECG) data ##
## -------------------------------------------------------------------------- ##

myGraph <- gfpop::graph(
  gfpop::Edge(0, 1, "down", penalty = 8e7, gap=0),
  gfpop::Edge(1, 2, "up", penalty = 0, gap=2000),
  gfpop::Edge(2, 3, "down", penalty = 0, gap=5000),
  gfpop::Edge(3, 4, "up", penalty = 0, gap=2000),
  gfpop::Edge(4, 5, "up", penalty = 0, gap=1000),
  gfpop::Edge(5, 6, "up", penalty = 0, gap=0),
  gfpop::Edge(6, 7, "down", penalty = 0, gap=0),
  gfpop::Edge(7, 8, "down", penalty = 0, gap=0),
  gfpop::Edge(8, 0, "up", penalty =0, gap=0),
                      all.null.edges = TRUE)

fit <- gfpop::gfpop(
  data = ECG$data$millivolts, mygraph = myGraph, type = "mean")

end.i <- fit$changepoints
start.i <- c(1, end.i[-length(end.i)]+1)
segments.dt <- with(ECG$data, data.table(
  timeStart=time[start.i],
  timeEnd=time[end.i],
  state=as.numeric(fit$states),
  mean=fit$parameters))

segments.dt[, letter := c("beforeQ", "Q", "R", "S", "S1", "S2", "peak", "afterPeak", "foo")[state+1] ]
mean.dt <- segments.dt[, data.table(
  time=as.numeric(rbind(timeStart-0.5, timeEnd+0.5)),
  mean=as.numeric(rbind(mean, mean)))]

m <- function(x){
  factor(x, c("previous", "changepoint"), c(
    "Previous model", "Proposed model"))
}
model.dt <- rbind(
  segments.dt[, data.table(
    model=m("changepoint"),
    time=ifelse(letter=="Q", timeEnd, (timeStart+timeEnd)/2),
    millivolts=mean,
    letter)],
  data.table(
    model=m("previous"),
    ECG$PanTompkins)
)[letter %in% c("Q", "R", "S")]

samples.per.second <- 250
truth.dt <- segments.dt[letter=="R", list(time=(timeStart+timeEnd)/2)]

gg <- ggplot()+
  geom_vline(aes(
    xintercept=time/samples.per.second),
    color="red",
    data=truth.dt)+
  geom_text(aes(
    x, y, hjust=hjust, label="True R"),
    color="red",
    size=3,
    data=data.table(
      x=208.5, y=6500, hjust=1, label="True R", model=m("changepoint")))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(model ~ .)+
  geom_line(aes(
    time/samples.per.second, millivolts),
    color="grey50",
    data=ECG$data)+
  geom_line(aes(
    time/samples.per.second, mean),
    data=data.table(model=m("changepoint"), mean.dt),
    color="blue")+
  geom_label(aes(
    time/samples.per.second, millivolts,
    label=letter),
    color="blue",
    size=3,
    label.padding=grid::unit(0.1, "lines"),
    alpha=0.6,
    data=model.dt)+
  coord_cartesian(xlim=c(52000, 52900)/samples.per.second, expand=FALSE)+
  xlab("Time (seconds)")+
  ylab("Electrocardiogram activity (mV)")

pdf( "figure13.pdf" )
print(gg)
dev.off()

