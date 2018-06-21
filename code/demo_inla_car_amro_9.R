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
  1 + log_hrs + std_yr +
  # cell unstructured and ICAR random effects
  f(grid_id1, model="besag", graph=g, scale.model=TRUE,
    hyper = list(prec = list(
      prior = "pc.prec",
      param = c(1, 0.01)))) +
  # cell unstructured and ICAR random effort slopes
  f(grid_id2, log_hrs, model="besag", graph=g, scale.model=TRUE,
    hyper = list(prec = list(
      prior = "pc.prec",
      param = c(1, 0.01)))) +
  # cell unstructured and ICAR random year slopes
  f(grid_id3, std_yr, model="besag", graph=g, scale.model=TRUE,
    hyper = list(prec = list(
      prior = "pc.prec",
      param = c(1, 0.01)))) +
  # random circle
  f(circle, model="iid",
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

# index and sort
modeling_data$grid_id2 <- modeling_data$grid_id1 <- modeling_data$grid_id
modeling_data$grid_id3 <- modeling_data$grid_id2
modeling_data <- arrange(modeling_data, grid_id, std_yr)



# run model and check output ---------------------------------------------------
out1<- inla(form1, family="nbinomial", data=modeling_data,
             # control.predictor=list(compute=T, link=1),
             control.compute=list(cpo=T, config=T),
             # control.inla=list(int.strategy='eb'),
             num.threads=3,
             verbose=T)

# view summaries
summary(out1, digits=3)
int <- out1$summary.fixed[1,1]
eft <- out1$summary.fixed[2,1]
yrs <- out1$summary.fixed[3,1]
g11 <- out1$summary.random$grid_id1$mean[1:n_cells]
g12 <- out1$summary.random$grid_id1$mean[(n_cells+1):(2*n_cells)]
alph <- exp(int + g11 + 0); hist(alph, breaks=60)
g21 <- out1$summary.random$grid_id2$mean[1:n_cells]
g22 <- out1$summary.random$grid_id2$mean[(n_cells+1):(2*n_cells)]
eps <- eft + g21 + 0; hist(eps)
g31 <- out1$summary.random$grid_id3$mean[1:n_cells]
g32 <- out1$summary.random$grid_id3$mean[(n_cells+1):(2*n_cells)]
tau <- yrs + g31 + 0; hist(tau)

# gof
sum(out1$cpo$failure, na.rm=T)
sum(log(out1$cpo$cpo), na.rm=T)
pit <- ggplot(data=data.frame(PIT=out1$cpo$pit), aes(x=PIT)) +
  geom_histogram(col="white") +
  xlab("Probability integral transform (PIT)") +
  ylab("Count")




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
epssamps <- ps1[grep("log_hrs", ps1$par_names), ][ , 1:posterior_ss]
row.names(epssamps) <- NULL
yrsamps <- ps1[grep("std_yr", ps1$par_names), ][ , 1:posterior_ss]
row.names(yrsamps) <- NULL
id1samps <- ps1[grep("id1", ps1$par_names), ][ , 1:posterior_ss][1:n_cells, ]
row.names(id1samps) <- NULL
id2samps <- ps1[grep("id2", ps1$par_names), ][ , 1:posterior_ss][1:n_cells, ]
row.names(id2samps) <- NULL
id3samps <- ps1[grep("id3", ps1$par_names), ][ , 1:posterior_ss][1:n_cells, ]
row.names(id3samps) <- NULL

# combine samples for cell alpha
al1 <- t(apply(as.matrix(id1samps), 1, function(x) x + as.numeric(intsamps)))
alpha <- exp(al1)
med_alpha <- apply(alpha, 1, FUN=median)
lcl_alpha <- apply(alpha, 1, FUN=quantile, probs=0.025)
ucl_alpha <- apply(alpha, 1, FUN=quantile, probs=0.975)
prec_alpha <- ucl_alpha - lcl_alpha
hist(med_alpha, 50)
hist(prec_alpha, 50)

# samples for epsilon
eps1 <- t(apply(as.matrix(id2samps), 1, function(x) x + as.numeric(epssamps)))
med_eps <- apply(eps1, 1, FUN=median)
lcl_eps <- apply(eps1, 1, FUN=quantile, probs=0.025)
ucl_eps <- apply(eps1, 1, FUN=quantile, probs=0.975)
prec_eps <- ucl_eps - lcl_eps
hist(as.matrix(eps1))
hist(prec_eps, 50)

# samples for tau
tau1 <- t(apply(as.matrix(id3samps), 1, function(x) x + as.numeric(yrsamps)))
tau1 <- (exp(tau1) - 1) * 100
med_tau <- apply(tau1, 1, FUN=median)
lcl_tau <- apply(tau1, 1, FUN=quantile, probs=0.025)
ucl_tau <- apply(tau1, 1, FUN=quantile, probs=0.975)
prec_tau <- ucl_tau - lcl_tau
hist(as.matrix(tau1))
hist(prec_tau, 50)

# collect posterior summaries into one dataframe
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
post_sum$bcr <- grid_key$bcr
post_sum$state <- grid_key$state
post_sum$country <- grid_key$country



# make cell level maps ---------------------------------------------------------
bcr_map <- readOGR(dsn=".", layer="simple_bcr")
bcr_sf <- as(bcr_map, "sf")
results_cells <- merge(modeling_cells, post_sum)
res_sf <- as(results_cells, "sf")

# filter out empty cells
hits1 <- which(!is.na(over(modeling_cells, sd2)["grid_id"])==T)
modeling_cells <- modeling_cells[hits1,]
plot(modeling_cells)
plot(sd2, add=T, pch=16, cex=0.5)
post_sum <- post_sum[hits1,]
summary(post_sum)

# map tau
tau1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=med_tau), col="gray40", size=0.3) +
  scale_fill_gradient2("Tau\n(% per year)", low = ("royalblue4"), mid = "white",
                         high = ("red4"), midpoint = 0, space = "Lab",
                         na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau2 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=prec_tau), col="gray40", size=0.3) +
  scale_fill_gradient2("Tau\ncredible\ninterval\nwidth\n(% per year)",
                       low = ("purple4"), mid = "white",
                       high = ("green4"), midpoint = 6, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau3 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=sig_tau), col="gray40", size=0.3) +
  scale_fill_gradient2("Significant tau\n(% per year)",
                       low = muted("royalblue4"), mid = "gray95",
                       high = muted("red4"), midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
#multiplot(tau1, tau3, tau2, cols=2)

# map epsilon
eps1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=med_eps), col="gray40", size=0.3) +
  scale_fill_gradient2("Epsilon", low = muted("purple4"), mid = "white",
                       high = muted("green4"), midpoint = 1, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
eps2 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=prec_eps), col="gray40", size=0.3) +
  scale_fill_gradient2("Epsilon\ninterval\nwidth", low = muted("purple4"),
                       mid = "white",
                       high = muted("orangered1"), midpoint = 0.5, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
eps3 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=sig_eps), col="gray40", size=0.3) +
  scale_fill_gradient2("Significant\nepsilon", low = muted("purple4"), mid = "white",
                       high = muted("green4"), midpoint = 1, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
#multiplot(eps1, eps3, cols=2)

# map alpha
alph1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=med_alpha*100), col="gray40", size=0.3) +
  scale_fill_gradient2("Alpha per\n100 hours", low = "tan4", mid = "white",
                       high = "green4", midpoint = 4.6, space = "Lab",
                       na.value = "grey40", guide = "colourbar", trans="log",
                       breaks=c(0.2, 5, 100, 2000)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
#multiplot(alph1, alph1, alph1, cols=2)

# print cell maps
pdf("tau_alpha.pdf", height=7.5, width=10.5)
multiplot(tau1, tau3, tau2, alph1, cols=2)
dev.off()

pdf("eps.pdf", height=3.75, width=10.5)
multiplot(eps1, eps3, cols=2)
dev.off()

pdf("eps_pit.pdf", height=3.75, width=10.5)
multiplot(eps1, pit, cols=2)
dev.off()



# aggregate results ------------------------------------------------------------
# summarise at bcr
new_bcr <- post_sum %>% group_by(bcr) %>%
  summarise(new_med_tau=median(med_tau), new_med_prec=median(prec_tau),
            new_min_prec=min(prec_tau), new_max_prec=max(prec_tau))
# bring in standard results
std_bcr <- read.csv("standard_trends.csv")
names(std_bcr) <- c("bcr", "std_med", "std_lcl", "std_ucl")
std_bcr$bcr <- gsub(pattern="BCR", replacement="", std_bcr$bcr)
std_bcr$std_prec <- std_bcr$std_ucl -std_bcr$std_lcl

# merge bcr results
compare_bcr_wide <- merge(bcr_sf, new_bcr, by.y="bcr", by.x="BCRNumber", all=F)
compare_bcr_wide <- merge(compare_bcr_wide, std_bcr,
                          by.y="bcr", by.x="BCRNumber", all=F)
summary(compare_bcr_wide)

# map side by side
names(compare_bcr_wide)
bcr_tau <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=compare_bcr_wide, aes(fill=new_med_tau), col="gray40", size=0.3) +
  scale_fill_gradient2("Tau\n(% per year)", low = "royalblue4", mid = "white",
                       high = "red4", midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar",
                       limits=c(-8, 15)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
bcr_trend <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=compare_bcr_wide, aes(fill=std_med), col="gray40", size=0.3) +
  scale_fill_gradient2("Standard\ntrend\n(% per year)", low = "royalblue4", mid = "white",
                       high = "red4", midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar",
                       limits=c(-8, 15)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# plot trends against each other
bcr_comp <- ggplot() +
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0, lty=2, col="gray70") +
  geom_vline(xintercept=0, lty=2, col="gray70") +
  xlim(c(-8, 15)) + ylim(c(-8, 15)) +
  ylab("Standard trend (% per year)") + xlab("Tau (% per year)") +
  geom_point(data=compare_bcr_wide, aes(x=new_med_tau, y=std_med),
             size=3, col="gray20")

# print cell maps
pdf("tau_bcr.pdf", height=7.5, width=10.5)
multiplot(bcr_tau, bcr_comp, bcr_trend, cols=2)
dev.off()

# correlation between trends
cor.test(x=compare_bcr_wide$new_med_tau, compare_bcr_wide$std_med, method="spearman")

# aggregate precision
compare_bcr_long <- as.data.frame(compare_bcr_wide) %>%
  select(1,5,6,7,11) %>%
  gather(key=Category, value=Precision, -BCRNumber)
compare_bcr_long$Category <- factor(compare_bcr_long$Category)
compare_bcr_long$Category <- factor(compare_bcr_long$Category,
                                    levels(compare_bcr_long$Category)[c(4,3,2,1)])

# histogram
svc_prec <- ggplot(data=compare_bcr_long) +
  geom_boxplot(aes(x=Category, y=Precision), col="gray20") +
  ylab("Credible interval width (% per year)") + xlab("Method") +
  scale_x_discrete(labels = c("std_prec" = "Standard",
                              "new_min_prec" = "SVC minimum",
                              "new_med_prec" = "SVC median",
                              "new_max_prec" = "SVC maximum"))

# std prec map
std_prec_map <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=compare_bcr_wide, aes(fill=std_prec), col="gray40", size=0.3) +
  scale_fill_gradient2("Standard\ninterval\nwidth\n(% per year)",
                       low = "green4", mid = "white",
                       high = "purple4", midpoint = 6, space = "Lab",
                       na.value = "grey40", guide = "colourbar",
                       limits=c(1, 12)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# new prec grid
new_prec_grid <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=prec_tau), col="gray40", size=0.3) +
  scale_fill_gradient2("SVC\ninterval\nwidth\n(% per year)", low = "green4", mid = "white",
                       high = "purple4", midpoint = 6, space = "Lab",
                       na.value = "grey40", guide = "colourbar",
                       limits=c(1, 20)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# print cell maps
pdf("prec_bcr.pdf", height=7.5, width=10.5)
multiplot(new_prec_grid, svc_prec, std_prec_map, cols=2)
dev.off()



# end analysis -----------------------------------------------------------------

write.csv(data.frame(sd2$lon, sd2$lat), "cbclatlon.csv")
