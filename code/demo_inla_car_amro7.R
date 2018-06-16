# cbc with car model -----------------------------------------------------------
# this is code for running an ICAR SVC model on CBC data.



# set up -----------------------------------------------------------------------
# install the required packages
library(RColorBrewer)
library(rgdal)
library(sp)
library(spdep)
library(sf)
library(concaveman)
library(raster)
library(tidyverse)
library(INLA)

# plot theme
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey80"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())

# theme
theme_map <- function(base_size = 9, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.background = element_blank(),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text.align=1,
          legend.text = element_text(size=9),
          legend.title=element_text(hjust=0, size=11),
          legend.justification=c(0, 0.5),
          plot.title = element_text(size=14, hjust = 0.7))
}

# multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# directories and options
wd1 <- "C:\\Users\\tmeehan\\Documents\\AudWork\\localCBCTrend\\spatialCBC\\AMRO"
wd2 <- "C:\\Users\\tmeehan\\Documents\\AudWork\\localCBCTrend\\spatialCBC"
options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)



# get and clean count data -----------------------------------------------------
setwd(wd1)
d1 <- read.csv("AMRO.csv")
d1 <- arrange(d1, strata, abbrev, count_yr)

# subset the data based on standard detection threshold criterion
d1 <- subset(d1, d1$DetectionsPerCircle > 0)
d1 <- subset(d1, d1$NonzeroCircles >= 5)
d1 <- subset(d1, d1$NumDetections / d1$NonzeroCircles >= 5)
d1 <- subset(d1, d1$TotalIndividuals >= 10)

# get rid of rediculous values
get_unlikely_vals <- function(dat, fld, stdevs=3, trans="log"){
  tfld <- which(names(dat)==fld)
  dat <- as.numeric(dat[,tfld])
  if(trans=="log") {dat <- log(dat)}
  if(trans=="log+1") {dat <- log(dat+1)}
  if(trans=="sqrt") {dat <- sqrt(dat)}
  scores <- scale(dat)
  idx <- which(scores < (-1 * stdevs) | scores > (1 * stdevs))
  return(idx)
}
d1 <- d1[-get_unlikely_vals(d1, "how_many", 3, "log+1"), ]
d1 <- d1[-get_unlikely_vals(d1, "hours", 3, "log"), ]

# sort and clean data
d1 <- arrange(d1, abbrev, bcr, count_yr)
d1[] <- lapply(d1, function(x) if(is.factor(x)) factor(x) else x)
d1[] <- lapply(d1, function(x) if(is.character(x)) as.factor(x) else x)

# get rid of junk
names(d1)
d2 <- data.frame(strata=d1$strata, circle=d1$abbrev, bcr=factor(d1$bcr),
                 circle_name=d1$circle_name, year=d1$count_yr,
                 count=d1$how_many,
                 hours=d1$hours, lon=d1$longitude, lat=d1$latitude,
                 area=d1$Area, detect_rat=d1$CircleDetectionRatio)

# map the points
crs4326 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs102008 <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
sd1 <- d2
coordinates(sd1) <- cbind(d2$lon, d2$lat)
proj4string(sd1) <- crs4326
sd1 <- spTransform(sd1, crs102008)

# swap in transformed coordinates
sd1$xcoord <- coordinates(sd1)[,1]
sd1$ycoord <- coordinates(sd1)[,2]

# get summaries
circles <- unique(as.character(sd1@data$circle))
n_circles <- length(circles) # 3234
strat <- unique(as.character(sd1@data$strata))
n_strat <- length(strat) # 156
bcrs <- unique(as.character(sd1@data$bcr))
n_bcrs <- length(bcrs) # 33
years <- unique(sd1@data$year)
n_years <- length(years) # 52
n_records <- n_circles * n_years # 168168

# add vars
d3 <- sd1@data
d3$yr_152 <- d3$year - min(d3$year) + 1
d3$std_yr <- d3$yr_152 - max(d3$yr_152)
d3$log_hrs <- log(d3$hours)
d3$obs <- row.names(d3)
d3$count <- as.integer(d3$count)
d3[] <- lapply(d3, function(x) if(is.character(x)) as.factor(x) else x)
d3[] <- lapply(d3, function(x) if(is.factor(x)) as.factor(x) else x)
str(d3)



# bring in na grid and match up with cbc circles -------------------------------
na_grid <- readOGR(dsn="../Shapes", "na_grid_w_centroids"); summary(na_grid)
proj4string(na_grid) <- crs102008

# make unique circle shapefile
d4 <- d3[ , c("strata", "circle", "circle_name", "lon", "lat",
             "xcoord", "ycoord")]
sd2 <- unique(d4[,])
row.names(sd2) <- NULL
coordinates(sd2) <- cbind(sd2$xcoord, sd2$ycoord)
proj4string(sd2) <- crs102008

# get concave hull of circle locations
hull_poly <- concaveman(sd2)

# grab grid cells covering circle loctions
circs <- st_transform(as(sd2, "sf"), 102008)
grids <- st_transform(as(na_grid, "sf"), 102008)
hulls <- st_transform(as(hull_poly, "sf"), 102008)
hulls <- st_buffer(hulls, 200000)
range1 <- as.integer(st_within(grids, hulls)) *
  (1:length(st_within(grids, hulls)))
range2 <- range1[!is.na(range1)]
range3 <- grids[range2, ]
range3 <- arrange(range3, as.numeric(id))
row.names(range3) <- NULL
range3$id <- as.numeric(row.names(range3))
cells_shape <- as(range3, "Spatial")
cells_shape <- cells_shape[order(cells_shape$id), ]
cells_shape$orig_id <- cells_shape$id
cells_shape$grid_id <- as.numeric(as.factor(cells_shape$id))
head(cells_shape@data)
plot(cells_shape, border="gray60", cex=0.5)
plot(hull_poly, add = TRUE, border="red")
plot(sd2, add=T, pch=16, col="blue", cex=0.3)
plot(cells_shape[c(1:150, 1150:nrow(cells_shape@data)), ], col="green", add=T)

# add grid id to unique circle layer
sd2$grid_id <- over(sd2, cells_shape)[,"grid_id"]
d5 <- merge(as.data.frame(sd2@data), as.data.frame(cells_shape@data),
            by="grid_id", by.y="grid_id", all.x=T, all.y=T, sort=F)
d5 <- dplyr::select(d5, circle, bcr_number, bcr_name, state=province, country,
                    grid_id, cell_lon=centroid_x, cell_lat=centroid_y)

# bring grid id to full data
d6 <- merge(d3, d5, by.x="circle", by.y="circle", all=T, sort=F)

# how many circles per cell
names(d6)
circles_per_cell <- d6 %>% dplyr::select(circle, grid_id, cell_lon, cell_lat) %>%
  unique() %>% arrange(grid_id) %>%
  group_by(grid_id, cell_lon, cell_lat) %>%
  summarise(number_circles=n())
ggplot(circles_per_cell, aes(number_circles)) +
  geom_histogram(breaks=c(0:20), fill="gray60", color="white") +
  xlab("Number of circles per grid cell") + ylab("Count")
ggplot(circles_per_cell, aes(x=cell_lon/1000, y=cell_lat/1000,
                             color=number_circles)) +
  geom_point(size=2, shape=15) +
  scale_color_distiller("Circles", palette="RdBu", na.value="gray80", trans="log",
                        breaks=c(0,1,2,4,8,16)) +
  xlab("Kilometers from center") + ylab("Kilometers from center")
summary(circles_per_cell)

# clean up and write merged data
modeling_data <- arrange((dplyr::select(d6, grid_id, circle, circle_name, strata,
                                        bcr=bcr_number, bcr_name=bcr_name,
                                        state, country, year, std_yr, count, hours,
                                        log_hrs, lon, lat, xcoord, ycoord)), grid_id, year)
modeling_data$obs <- 1:nrow(modeling_data)
modeling_data$log_hrs[is.na(modeling_data$log_hrs)] <- 0
modeling_data$std_yr[is.na(modeling_data$std_yr)] <- 0
write.csv(modeling_data, "modeling_data.csv", row.names=F)
modeling_cells <- cells_shape
modeling_cells@data <- select(modeling_cells@data, grid_id, cell_x=centroid_x,
                              cell_y=centroid_y)
modeling_cells <- modeling_cells[order(modeling_cells$grid_id), ]
row.names(modeling_cells@data) <- modeling_cells@data$grid_id
head(modeling_cells@data); n_cells <- nrow(modeling_cells@data)
writeOGR(modeling_cells, dsn=getwd(), layer="modeling_cells",
         driver="ESRI Shapefile")

# make grid key for later
grid_key <- unique(modeling_data[, c("grid_id", "bcr", "bcr_name", "state",
                                     "country")])
row.names(grid_key) <- NULL



# make neighbors ---------------------------------------------------------------
nb1 <- poly2nb(modeling_cells, row.names=modeling_cells$grid_id); nb1
is.symmetric.nb(nb1, verbose = FALSE, force = TRUE)
nb2INLA("nb1.graph", nb1)
nb1.adj <- paste(getwd(),"/nb1.graph", sep="")
g <- inla.read.graph("nb1.graph");  # plot(g)
head(g$nbs)
plot(nb1, coordinates(modeling_cells), col="red", cex=0.5)



# specify model and prep data --------------------------------------------------
form1 <- count ~
  # grand mean
  1 +
  # cell unstructured and ICAR random effects
  f(grid_id1, model="bym2", constr=T, graph=g, scale.model=TRUE,
    hyper = list(prec = list(
      prior = "pc.prec",
      param = c(1, 0.01)))) +
  # cell unstructured and ICAR random effort slopes
  f(grid_id2, log_hrs, model="besag", graph=g, constr=F, scale.model=TRUE,
    hyper = list(prec = list(
      prior = "pc.prec",
      param = c(1, 0.01)))) +
  # cell unstructured and ICAR random year slopes
  f(grid_id3, std_yr, model="besag", graph=g, constr=F, scale.model=TRUE,
    hyper = list(prec = list(
      prior = "pc.prec",
      param = c(1, 0.01)))) +
  # overdispersion effect
  f(obs, model="iid", constr=T)

# index and sort
modeling_data$grid_id2 <- modeling_data$grid_id1 <- modeling_data$grid_id
modeling_data$grid_id3 <- modeling_data$grid_id2
modeling_data <- arrange(modeling_data, grid_id, std_yr)



# run model and check output ---------------------------------------------------
out1<- inla(form1, family="poisson", data=modeling_data,
             control.predictor=list(compute=T, link=1),
             control.compute=list(dic=F, cpo=T, waic=T, config=T),
             control.inla=list(int.strategy='eb'),
             num.threads=2,
             verbose=T)

# view summaries
summary(out1, digits=3)
int <- out1$summary.fixed$mean
g11 <- out1$summary.random$grid_id1$mean[1:n_cells]
g12 <- out1$summary.random$grid_id1$mean[(n_cells+1):(2*n_cells)]
alph <- exp(int + g11 + 0); hist(alph, breaks=60)
g21 <- out1$summary.random$grid_id2$mean[1:n_cells]
g22 <- out1$summary.random$grid_id2$mean[(n_cells+1):(2*n_cells)]
eps <- g21 + 0; hist(eps)
g31 <- out1$summary.random$grid_id3$mean[1:n_cells]
g32 <- out1$summary.random$grid_id3$mean[(n_cells+1):(2*n_cells)]
tau <- g31 + 0; hist(tau)

# gof
sum(out1$cpo$failure, na.rm=T)
hist(out1$cpo$pit)


# take samples from the marginal posteriors ------------------------------------
# get sample
posterior_ss <- 500
samp <- inla.posterior.sample(posterior_ss, out1, num.threads=3)
par_names <- as.character(attr(samp[[1]]$latent, "dimnames")[[1]])
ps1 <- as.data.frame(sapply(samp, function(x) x$latent))
ps1$par_names <- par_names
rm(samp)

# break it up for cell analysis
intsamps <- ps1[grep("(Intercept)", ps1$par_names), ][ , 1:posterior_ss]
row.names(intsamps) <- NULL
id1samps <- ps1[grep("id1", ps1$par_names), ][ , 1:posterior_ss][1:n_cells, ]
row.names(id1samps) <- NULL
id2samps <- ps1[grep("id2", ps1$par_names), ][ , 1:posterior_ss][1:n_cells, ]
row.names(id2samps) <- NULL
id3samps <- ps1[grep("id3", ps1$par_names), ][ , 1:posterior_ss][1:n_cells, ]
row.names(id3samps) <- NULL

# combine samples for cell alpha
al1 <- t(apply(as.matrix(id1samps), 1, function(x) x +
                 as.numeric(intsamps)))
alpha <- exp(al1)
med_alpha <- apply(alpha, 1, FUN=median)
lcl_alpha <- apply(alpha, 1, FUN=quantile, probs=0.025)
ucl_alpha <- apply(alpha, 1, FUN=quantile, probs=0.975)
prec_alpha <- ucl_alpha - lcl_alpha
hist(med_alpha, 50)

# samples for epsilon
hist(as.matrix(id2samps))
med_eps <- apply(id2samps, 1, FUN=median)
lcl_eps <- apply(id2samps, 1, FUN=quantile, probs=0.025)
ucl_eps <- apply(id2samps, 1, FUN=quantile, probs=0.975)
prec_eps <- ucl_eps - lcl_eps

# samples for tau
hist(as.matrix(id3samps))
tau1 <- (exp(id3samps) - 1) * 100
med_tau <- apply(tau1, 1, FUN=median)
lcl_tau <- apply(tau1, 1, FUN=quantile, probs=0.025)
ucl_tau <- apply(tau1, 1, FUN=quantile, probs=0.975)
prec_tau <- ucl_tau - lcl_tau

# collect posterior summaries
post_sum <- data.frame(grid_id=1:length(med_alpha),
                       med_alpha, lcl_alpha, ucl_alpha, prec_alpha,
                       med_eps, lcl_eps, ucl_eps, prec_eps, sig_eps=NA,
                       med_tau, lcl_tau, ucl_tau, prec_tau, sig_tau=NA)
post_sum$sig_eps <- ifelse((post_sum$lcl_eps < 1 & post_sum$ucl_eps > 1),
                           post_sum$sig_eps <- NA,
                           post_sum$sig_eps <- post_sum$med_eps)
post_sum$sig_tau <- ifelse((post_sum$lcl_tau < 1 & post_sum$ucl_tau > 0),
                           post_sum$sig_tau <- NA,
                           post_sum$sig_tau <- post_sum$med_tau)



# make cell level maps ---------------------------------------------------------
bcr_map <- readOGR(dsn="../Shapes", layer="simple_bcr")
bcr <- as(bcr_map, "sf")
results_cells <- merge(modeling_cells, post_sum)
res <- as(results_cells, "sf")

# map tau
tau1 <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=res, aes(fill=med_tau), col="gray40") +
  scale_fill_distiller("Tau", palette="RdYlBu", direction=1,
                       breaks=seq(-10, 10, 5), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau2 <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=res, aes(fill=prec_tau), col="gray40") +
  scale_fill_distiller("Interval\nwidth", palette="YlGnBu", direction=1,
                       breaks=seq(0, 15, 5), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau3 <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=res, aes(fill=sig_tau), col="gray40") +
  scale_fill_distiller("Significant\ntau", palette="RdYlBu", direction=1,
                       breaks=seq(-10, 10, 5), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
multiplot(tau1, tau3, tau2, cols=2)

# map epsilon
eps1 <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=res, aes(fill=med_eps), col="gray40") +
  scale_fill_distiller("Epsilon", palette="PRGn", direction=1,
                       breaks=seq(0, 2.4, 0.8), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
eps2 <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=res, aes(fill=prec_eps), col="gray40") +
  scale_fill_distiller("Interval\nwidth", palette="YlGnBu", direction=1,
                       breaks=seq(0, 3, 1), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
eps3 <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=res, aes(fill=sig_eps), col="gray40") +
  scale_fill_distiller("Significant\nepsilon", palette="PRGn", direction=1,
                       breaks=seq(0, 2.4, 0.8), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
multiplot(eps1, eps3, eps2, cols=2)

# map alpha
alph1 <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=res, aes(fill=med_alpha*100), col="gray40") +
  scale_fill_distiller("Alpha per\n100 hours", palette="BrBG",
                       direction=1,trans="log",
                       breaks=c(0.2, 5, 100, 2000), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
multiplot(alph1, alph1, alph1, cols=2)



# scale up to BCR and compare --------------------------------------------------
# merge aggregation info
n_keys <- length(grid_key)
intsamps_ag <- as.data.frame(cbind(grid_key, intsamps))
id1samps_ag <- as.data.frame(cbind(grid_key, id1samps))
id2samps_ag <- as.data.frame(cbind(grid_key, id2samps))
id3samps_ag <- as.data.frame(cbind(grid_key, id3samps))

# # filter for observed cells
# hits1 <- which(!is.na(over(modeling_cells, sd2)["grid_id"])==T)
# intsamps_ag <- intsamps_ag[hits1,]
# id1samps_ag <- id1samps_ag[hits1,]
# id2samps_ag <- id2samps_ag[hits1,]
# id3samps_ag <- id3samps_ag[hits1,]

# process alpha samples
bcrs <- unique(id1samps_ag$bcr)
bcr_par_sum <- NULL
for(i in 1:length(bcrs)){
  grab1 <- which(id1samps_ag$bcr==as.numeric(bcrs[i]))
  l1 <- as.matrix(id1samps_ag[grab1, c((n_keys+1):(n_keys+posterior_ss))])
  al1 <- t(apply(l1, 1, function(x) x + as.numeric(intsamps)))
  al2 <- c(exp(al1))
  med_alpha <- median(al2, na.rm=T)
  lcl_alpha <- quantile(al2, probs=0.025)
  ucl_alpha <- quantile(al2, probs=0.975)
  prec_alpha <- ucl_alpha - lcl_alpha
  c1 <- data.frame(bcr=bcrs[i], par="alpha", med=med_alpha, lcl=lcl_alpha,
                   ucl=ucl_alpha, prec=prec_alpha)
  row.names(c1) <- NULL
  bcr_par_sum <- rbind(bcr_par_sum, c1)

}
for(i in 1:length(bcrs)){
  grab1 <- which(id2samps_ag$bcr==as.numeric(bcrs[i]))
  l1 <- as.matrix(id2samps_ag[grab1, c((n_keys+1):(n_keys+posterior_ss))])
  al2 <- c(l1)
  med_eps <- median(al2, na.rm=T)
  lcl_eps <- quantile(al2, probs=0.025)
  ucl_eps <- quantile(al2, probs=0.975)
  prec_eps <- ucl_eps - lcl_eps
  c1 <- data.frame(bcr=bcrs[i], par="eps", med=med_eps, lcl=lcl_eps,
                   ucl=ucl_eps, prec=prec_eps)
  row.names(c1) <- NULL
  bcr_par_sum <- rbind(bcr_par_sum, c1)
}
for(i in 1:length(bcrs)){
  grab1 <- which(id3samps_ag$bcr==as.numeric(bcrs[i]))
  l1 <- as.matrix(id3samps_ag[grab1, c((n_keys+1):(n_keys+posterior_ss))])
  al1 <- (exp(l1) - 1) * 100
  al2 <- al1
  med_tau <- median(al2, na.rm=T)
  lcl_tau <- quantile(al2, probs=0.025)
  ucl_tau <- quantile(al2, probs=0.975)
  prec_tau <- ucl_tau - lcl_tau
  c1 <- data.frame(bcr=bcrs[i], par="tau", med=med_tau, lcl=lcl_tau,
                   ucl=ucl_tau, prec=prec_tau)
  row.names(c1) <- NULL
  bcr_par_sum <- rbind(bcr_par_sum, c1)
}
bcr_tau <- filter(bcr_par_sum, par=="tau")

# merge with bcr map
bcr_map <- merge(bcr_map, bcr_tau, by.x="BCRNumber", by.y="bcr", all=T)
bcr_sf <- as(bcr_map, "sf")
ggplot() +
  geom_sf(data=bcr_sf, aes(fill=med), col="gray40") +
  scale_fill_distiller("Tau", palette="RdYlBu", direction=1,
                       breaks=seq(-10, 10, 5), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
ggplot() +
  geom_sf(data=bcr_sf, aes(fill=prec), col="gray40") +
  scale_fill_distiller("Interval", palette="RdYlBu", direction=1,
                       breaks=seq(0, 15, 5), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# bring in standard results
std_bcr <- read.csv("standard_trends.csv")
names(std_bcr) <- c("BCRNumber", "mean_trend", "lcl", "ucl")
std_bcr$BCRNumber <- gsub(pattern="BCR", replacement="", std_bcr$BCRNumber)
bcr_com <- merge(bcr_sf, std_bcr, by="BCRNumber")
bcr_com$prec.y <- bcr_com$ucl.y - bcr_com$lcl.y

# map side by side
bcr_tau <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=bcr_com, aes(fill=med), col="gray40") +
  scale_fill_distiller("Tau", palette="RdYlBu", direction=1,
                       breaks=seq(-8, 16, 4), limits=c(-8, 16), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
bcr_trend <- ggplot() +
  geom_sf(data=bcr, fill="gray40", col="gray40") +
  geom_sf(data=bcr_com, aes(fill=mean_trend), col="gray40") +
  scale_fill_distiller("Trend", palette="RdYlBu", direction=1,
                       breaks=seq(-8, 16, 4), limits=c(-8, 16), na.value="gray40") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# plot trends against each other
bcr_comp <- ggplot() +
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0, lty=2, col="gray70") +
  geom_vline(xintercept=0, lty=2, col="gray70") +
  xlim(c(-8, 15)) + ylim(c(-8, 15)) +
  xlab("BCR trend, standard method") + ylab("BCR trend, SVC model") +
  geom_point(data=bcr_com, aes(x=mean_trend, y=med), size=3, col="gray20")
multiplot(bcr_trend, bcr_comp, bcr_tau, cols=2)

# comp stats
cor.test(x=bcr_com$med, bcr_com$mean_trend, method="spearman")



# end analysis


