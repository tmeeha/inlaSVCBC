# cbc with car model -----------------------------------------------------------
# this is code for running an ICAR SVC model on CBC data.
# ------------------------------------------------------------------------------



# set up -----------------------------------------------------------------------
# install the required packages
library(RColorBrewer)
library(tidyverse)
library(rgdal)
library(sp)
library(spdep)
library(sf)
library(concaveman)
library(scales)
library(INLA)
library(brinla)
library(inlabru)

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

# map theme
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

# directories and options
wd1 <- "C:\\Users\\tmeehan\\Documents\\Github\\inlaSVCBC\\code"
wd2 <- "C:\\Users\\tmeehan\\Documents\\Github\\inlaSVCBC\\code"
options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)
# ------------------------------------------------------------------------------



# get count data and grids -----------------------------------------------------
setwd(wd1)
modeling_data <- read.csv("modeling_data.csv"); head(modeling_data)
cbc_amro_circles <- readOGR(dsn=".", layer="cbc_amro_circles")
head(cbc_amro_circles@data)
cbc_amro_grid <- readOGR(dsn=".", layer="cbc_amro_grid")
names(cbc_amro_grid)[9] <- "grid_id"; head(cbc_amro_grid@data)
cbc_na_grid <- readOGR(dsn=".", layer="cbc_na_grid")
cbc_na_grid$na_id <- cbc_na_grid$id; head(cbc_na_grid@data)
plot(cbc_na_grid)
plot(cbc_amro_grid, col="red", add=T)
grid1 <- as(cbc_amro_grid, "sf")
plot(cbc_amro_circles, col="blue", add=T)
grid_key <- unique(modeling_data[, c("grid_id", "bcr", "bcr_name", "state",
                                     "country")])
row.names(grid_key) <- NULL
# ------------------------------------------------------------------------------



# summarize circles per grid ---------------------------------------------------
# how many circles per all 1316 cells
names(modeling_data)
circles_per_cell <- modeling_data %>%
  dplyr::select(circle, grid_id, lon, lat) %>%
  arrange(grid_id) %>% unique() %>%
  mutate(occupied=ifelse(is.na(lon), 0, 1)) %>%
  group_by(grid_id) %>%
  summarise(number_circles=sum(occupied))
ggplot(circles_per_cell, aes(number_circles)) +
  geom_histogram(breaks=c(0:20), fill="gray60", color="white") +
  xlab("Number of circles per grid cell") + ylab("Count")
summary(circles_per_cell)

# how many circles per 880 occupied cells
circles_per_cell <- modeling_data %>%
  dplyr::select(circle, grid_id, lon, lat) %>%
  arrange(grid_id) %>% unique() %>%
  mutate(occupied=ifelse(is.na(lon), 0, 1)) %>%
  filter(occupied==1) %>%
  group_by(grid_id) %>%
  summarise(number_circles=sum(occupied))
ggplot(circles_per_cell, aes(number_circles)) +
  geom_histogram(breaks=c(0:20), fill="gray60", color="white") +
  xlab("Number of circles per grid cell") + ylab("Count")
summary(circles_per_cell)
# ------------------------------------------------------------------------------



# make neighbors ---------------------------------------------------------------
nb1 <- poly2nb(cbc_amro_grid, row.names=cbc_amro_grid$grid_id); nb1
is.symmetric.nb(nb1, verbose = FALSE, force = TRUE)
nb2INLA("nb1.graph", nb1)
nb1.adj <- paste(getwd(),"/nb1.graph", sep="")
g1 <- inla.read.graph("nb1.graph")
plot(nb1, coordinates(cbc_amro_grid), col="red", cex=0.5)

# index and sort
modeling_data$eps_i <- modeling_data$alpha_i <- modeling_data$grid_id
modeling_data$tau_i <- modeling_data$eps_i
modeling_data$kappa_k <- as.integer(factor(modeling_data$circle))
modeling_data <- arrange(modeling_data, grid_id, std_yr)
n_circs <- max(modeling_data$kappa_k, na.rm=T)
n_cells <- max(modeling_data$alpha_i, na.rm=T)
# ------------------------------------------------------------------------------



# make and run model -----------------------------------------------------------
# model
form1 <- count ~ -1 + # remove grand mean
  # cell ICAR random intercepts
  f(alpha_i, model="besag", graph=g1, constr=FALSE, scale.model=TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # cell ICAR random effort slopes
  f(eps_i, log_hrs, model="besag", graph=g1, constr=FALSE, scale.model=TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # cell ICAR random year slopes
  f(tau_i, std_yr, model="besag", graph=g1, constr=FALSE, scale.model=TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  # random circle intercepts
  f(kappa_k, model="iid", constr=TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

# get results
out1<- inla(form1, family="nbinomial", data=modeling_data,
            control.compute=list(cpo=T, config=T),
            control.inla=list(strategy="adaptive",
                              int.strategy="auto"),
            num.threads=3,
            verbose=T)

# view summaries
summary(out1, digits=3)
bri.hyperpar.summary(out1)

# remove cells with no circles
cells_with_counts <- unique(modeling_data$grid_id[which(
  !is.na(modeling_data$hours))])

# get alpha summaries
alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
alph_iw <- alph_ul - alph_ll
hist(alph); summary(alph)
hist(alph_ll); summary(alph_ll)
hist(alph_ul); summary(alph_ul)

# get epsilon summaries
eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
eps_iw <- eps_ul - eps_ll
hist(eps); summary(eps); round(sum(eps<1)/length(eps), 2)
hist(eps_ll); summary(eps_ll)
hist(eps_ul); summary(eps_ul)
round(sum(eps_ll<=1 & eps_ul>=1)/length(eps_ll), 2)
round(sum(eps_ll<=0 & eps_ul>=0)/length(eps_ll), 2)
cor.test(eps, alph, method="spearman")

# get tau summaries
tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts])
        - 1) * 100
tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts])
           - 1) * 100
tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts])
           - 1) * 100
tau_iw <- tau_ul - tau_ll
hist(tau); summary(tau); round(sum(tau>=0)/length(tau), 2)
hist(tau_ll); summary(tau_ll)
hist(tau_ul); summary(tau_ul)
hist(tau_iw); summary(tau_iw)
round(sum(tau_ll<=0 & tau_ul<=0)/length(tau_ll), 2)
round(sum(tau_ll>=0 & tau_ul>=0)/length(tau_ll), 2)
cor(cbind(alph, tau), method="spearman")
ppcor::pcor(cbind(eps, alph, tau), method="spearman")

# gof
sum(out1$cpo$failure, na.rm=T)
-2 * sum(log(out1$cpo$cpo[out1$cpo$failure==0]), na.rm=T)
pit1 <- data.frame(PIT=out1$cpo$pit) %>%
  filter(out1$cpo$pit<0.99 & out1$cpo$failure!=1 & out1$cpo$pit>0.01)
pit2 <- ggplot(data=pit1, aes(x=PIT)) +
  geom_histogram(col="white") +
  xlab("Probability integral transform (PIT)") +
  ylab("Count"); pit2; summary(pit1$PIT)
# ------------------------------------------------------------------------------



# # specify model with year cell effects and prep data -------------------------
# modeling_data$gamma_ij <- paste0(modeling_data$grid_id, "-",
#   modeling_data$year)
#
# form2 <- count ~ -1 + # remove grand mean
#   # cell ICAR random intercepts
#   f(alpha_i, model="besag", graph=g1, constr=FALSE, scale.model=TRUE,
#     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#   # cell ICAR random effort slopes
#   f(eps_i, log_hrs, model="besag", graph=g1, constr=FALSE, scale.model=TRUE,
#     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#   # cell ICAR random year slopes
#   f(tau_i, std_yr, model="besag", graph=g1, constr=FALSE, scale.model=TRUE,
#     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#   # random circle intercepts
#   f(kappa_k, model="iid", constr=TRUE,
#     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#   # cell-year effect
#   f(gamma_ij, model="iid", constr=TRUE,
#     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
#
# out2<- inla(form2, family="nbinomial", data=modeling_data,
#             control.compute=list(cpo=T, config=T),
#             control.inla=list(strategy="adaptive",
#                               int.strategy='auto'),
#             num.threads=3,
#             verbose=T)
#
# # view summaries
# summary(out2, digits=3)
# bri.hyperpar.summary(out2)
# sum(out1$cpo$failure, na.rm=T)
# sum(out2$cpo$failure, na.rm=T)
# -2 * mean(log(out1$cpo$cpo[out1$cpo$failure==0]), na.rm=T)
# -2 * mean(log(out2$cpo$cpo[out2$cpo$failure==0]), na.rm=T)
# pit1 <- data.frame(PIT=out2$cpo$pit) %>%
#   filter(out2$cpo$pit<0.99 & out2$cpo$failure!=1 & out2$cpo$pit>0.01)
# pit2 <- ggplot(data=pit1, aes(x=PIT)) +
#   geom_histogram(col="white") +
#   xlab("Probability integral transform (PIT)") +
#   ylab("Count"); pit2; summary(pit1$PIT)
#
# # time series plots per cell
# cell_ts <- function(cell1=801){
#   d0 <- out2$summary.random$alpha_i$`0.5quant`[cell1]
#   d1 <- out2$summary.random$tau_i$`0.5quant`[cell1]
#   d2 <- data.frame(
#   year=as.numeric(gsub(paste0(cell1,"-"), "",
#                        grep(paste0("\\b",cell1,"-"),
#                             out2$summary.random$gamma_ij$ID,
#                                  value=T)))-117,
#   gamma_ij=
#     out2$summary.random$gamma_ij$`0.5quant`[grep(
#       paste0("\\b",cell1,"-"), out2$summary.random$gamma_ij$ID)]) %>%
#     arrange(year)
#   d2$x0 <- d0
#   d2$x1 <- d2$year*d1
#   d2$abund <- exp(d2$x0 + d2$x1 + d2$gamma_ij)
#   d2$trend <- exp(d2$x0 + d2$x1)
#   ggplot(d2, aes(x=year+117+ 1900)) + geom_line(aes(y=trend)) +
#   geom_point(aes(y=abund)) + annotate("text", x=1975, y=14,
#                                       label=paste0("Grid cell number ", cell1))
# }
# cell2 <- c(798, 799, 800, 801)
# for(i in cell2){
#   print(cell_ts(i))
# }
#
#
# # ----------------------------------------------------------------------------



# play with samples ------------------------------------------------------------
posterior_ss <- 10
samp1 <- inla.posterior.sample(posterior_ss, out1, num.threads=3)
par_names <- as.character(attr(samp1[[1]]$latent, "dimnames")[[1]])
post1 <- as.data.frame(sapply(samp1, function(x) x$latent))
post1$par_names <- par_names

# tau samples
tau_samps1 <- post1[grep("tau_i", post1$par_names), ]
row.names(tau_samps1) <- NULL
tau_samps1 <- tau_samps1[cells_with_counts, 1:posterior_ss]
tau_samps1 <- (exp(tau_samps1) - 1) * 100
tau_samps2 <- cbind(grid_key[cells_with_counts,], tau_samps1)
row.names(tau_samps2) <- NULL
val_names <- grep("V", names(tau_samps2))

# tau_bcr
tau_bcr <- tau_samps2 %>%
  dplyr::select(bcr, val_names) %>%
  mutate(bcr=factor(bcr)) %>%
  gather(key=key, val=val, -bcr) %>%
  dplyr::select(-key) %>%
  group_by(bcr) %>%
  summarise(med_tau=median(val), lcl_tau=quantile(val, probs=0.025),
            ucl_tau=quantile(val, probs=0.975), iw_tau=ucl_tau-lcl_tau,
            n=n()/posterior_ss); View(tau_bcr)

# tau_total
tau_tot <- tau_samps2 %>%
  dplyr::select(val_names) %>%
  as.matrix() %>%
  as.numeric() %>%
  quantile(probs=c(0.5, 0.025, 0.975)); tau_tot
# ------------------------------------------------------------------------------



# collect posterior summaries into one dataframe -------------------------------
post_sum <- data.frame(grid_id=cells_with_counts,
                       alph, alph_ll, alph_ul, alph_iw,
                       eps, eps_ll, eps_ul, eps_iw, eps_sig=NA,
                       tau, tau_ll, tau_ul, tau_iw, tau_sig=NA)
post_sum$eps_sig <- ifelse((post_sum$eps_ll < 1 & post_sum$eps_ul > 1),
                           post_sum$eps_sig <- NA,
                           post_sum$eps_sig <- post_sum$eps)
post_sum$tau_sig <- ifelse((post_sum$tau_ll < 1 & post_sum$tau_ul > 1),
                           post_sum$tau_sig <- NA,
                           post_sum$tau_sig <- post_sum$tau)
post_sum$bcr <- grid_key$bcr[cells_with_counts]
post_sum$state <- grid_key$state[cells_with_counts]
post_sum$country <- grid_key$country[cells_with_counts]
summary(post_sum)
# ------------------------------------------------------------------------------



# make cell level maps ---------------------------------------------------------
bcr_map <- readOGR(dsn=".", layer="simple_bcr")
bcr_sf <- as(bcr_map, "sf")
results_cells <- merge(cbc_amro_grid, post_sum, by="grid_id", all=F)
res_sf <- as(results_cells, "sf")
plot(res_sf["tau"])
cbc_amro_grid <- cbc_amro_grid[cells_with_counts, ]
plot(cbc_amro_grid)
plot(cbc_amro_circles, add=T, pch=16, cex=0.5)

# grid
circles1 <- as(cbc_amro_circles, "sf")
grid1 <- merge(grid1, circles_per_cell, all=T)
grid1$number_circles[is.na(grid1$number_circles)] <- 0
grid_p1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=grid1, aes(fill=number_circles), col="gray40", size=0.3) +
  scale_fill_gradient("A. Circles\nper cell", low = ("white"),
                       high = ("red4"), space = "Lab", trans="log1p",
                       na.value = "grey40", guide = "colourbar",
                      breaks=c(0,2,7,20)) +
  theme_map() +
  theme(panel.grid.major=element_line(colour="transparent"))

# map tau
tau_p1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau), col="gray40", size=0.3) +
  scale_fill_gradient2("D. Tau\n(% per year)", low = ("royalblue4"),
                       mid = "white",
                       high = ("red4"), midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau_p2 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_iw), col="gray40", size=0.3) +
  scale_fill_gradient2("F. Tau\ncredible\ninterval\nwidth\n(% per year)",
                       low = ("purple4"), mid = "white",
                       high = ("green4"), midpoint = 6, space = "Lab",
                       na.value = "grey40", guide = "colourbar",
                       breaks=c(3,6,9,12)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau_p3 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_sig), col="gray40", size=0.3) +
  scale_fill_gradient2("E. Significant\ntau (% per year)",
                       low = muted("royalblue4"), mid = "gray95",
                       high = muted("red4"), midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
# multiplot(tau_p1)
# multiplot(tau_p1, tau_p3, tau_p2, cols=2)

# map epsilon
eps_p1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=eps), col="gray40", size=0.3) +
  scale_fill_gradient2("C. Epsilon", low = muted("purple4"), mid = "white",
                       high = muted("green4"), midpoint = 1, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
eps_p2 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=eps_iw), col="gray40", size=0.3) +
  scale_fill_gradient2("Epsilon\ninterval\nwidth", low = muted("purple4"),
                       mid = "white",
                       high = muted("orangered1"), midpoint = 0.5,
                       space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
eps_p3 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=eps_sig), col="gray40", size=0.3) +
  scale_fill_gradient2("Significant\nepsilon", low = muted("purple4"),
                       mid = "white",
                       high = muted("green4"), midpoint = 1, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
# multiplot(eps_p1)
# multiplot(eps_p1, eps_p3, eps_p2, cols=2)

# map alpha
alph_p1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=alph), col="gray40", size=0.3) +
  scale_fill_gradient2("B. Alpha", low = "tan4", mid = "white",
                       high = "green4", midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar", trans="log",
                       breaks=c(0.01, 0.1, 1, 10)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
# multiplot(alph_p1)

# print cell maps
pdf("tau_alpha.pdf", height=11.25, width=10.5)
multiplot(grid_p1, eps_p1, tau_p3, alph_p1, tau_p1, tau_p2,cols=2)
dev.off()

pdf("eps.pdf", height=3.75, width=10.5)
multiplot(eps_p1, eps_p3, cols=2)
dev.off()
# ------------------------------------------------------------------------------



# aggregate results ------------------------------------------------------------
# summarise at bcr
new_bcr <- post_sum %>% group_by(bcr) %>%
  summarise(new_med_tau=median(tau), new_med_prec=median(tau_iw),
            new_min_prec=min(tau_iw), new_max_prec=max(tau_iw))
# bring in standard results
std_bcr <- read.csv("standard_trends.csv")
names(std_bcr) <- c("bcr", "std_med", "std_lcl", "std_ucl")
std_bcr$bcr <- gsub(pattern="BCR", replacement="", std_bcr$bcr)
std_bcr$std_prec <- std_bcr$std_ucl - std_bcr$std_lcl

# merge bcr results
compare_bcr_wide <- merge(bcr_sf, new_bcr, by.y="bcr", by.x="BCRNumber",
                          all=F)
compare_bcr_wide <- merge(compare_bcr_wide, std_bcr,
                          by.y="bcr", by.x="BCRNumber", all=F)
summary(compare_bcr_wide)
sd(compare_bcr_wide$new_med_prec)
sd(compare_bcr_wide$std_prec)

# map side by side
names(compare_bcr_wide)
bcr_tau <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=compare_bcr_wide, aes(fill=new_med_tau), col="gray40",
          size=0.3) +
  scale_fill_gradient2("Tau\n(% per year)", low = "royalblue4",
                       mid = "white",
                       high = "red4", midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar",
                       limits=c(-8, 15)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
bcr_trend <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=compare_bcr_wide, aes(fill=std_med), col="gray40",
          size=0.3) +
  scale_fill_gradient2("Standard\ntrend\n(% per year)", low = "royalblue4",
                       mid = "white",
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

# print bcr maps
pdf("tau_bcr.pdf", height=7.5, width=10.5)
multiplot(bcr_tau, bcr_comp, bcr_trend, cols=2)
dev.off()

# correlation between trends
cor.test(x=compare_bcr_wide$new_med_tau, compare_bcr_wide$std_med,
         method="spearman")

# aggregate precision
compare_bcr_long <- as.data.frame(compare_bcr_wide) %>%
  dplyr::select(1,5,6,7,11) %>%
  gather(key=Category, value=Precision, -BCRNumber)
compare_bcr_long$Category <- factor(compare_bcr_long$Category)
compare_bcr_long$Category <- factor(compare_bcr_long$Category,
                                    levels(compare_bcr_long$Category)
                                    [c(4,3,2,1)])
compare_bcr_long %>% group_by(Category) %>% summarise(mean(Precision))
coin::independence_test(Precision~Category,
             data=subset(compare_bcr_long,
                         Category!="new_min_prec" & Category!="new_max_prec"))

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
                       high = "green4", mid = "white",
                       low = "purple4", midpoint = 6, space = "Lab",
                       na.value = "grey40", guide = "colourbar",
                       limits=c(1, 20)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# new prec grid
new_prec_grid <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_iw), col="gray40", size=0.3) +
  scale_fill_gradient2("SVC\ninterval\nwidth\n(% per year)", high = "green4",
                       mid = "white",
                       low = "purple4", midpoint = 6, space = "Lab",
                       na.value = "grey40", guide = "colourbar",
                       limits=c(1, 20)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))

# print cell maps
pdf("prec_bcr.pdf", height=7.5, width=10.5)
multiplot(new_prec_grid, svc_prec, std_prec_map, cols=2)
dev.off()
# ------------------------------------------------------------------------------



# import trend covariates ------------------------------------------------------
trend_covs <- readOGR("cbc_amro_grid_covs.shp")
names(trend_covs)
trend_covs <- spTransform(trend_covs, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
trend_covs$lon <- coordinates(trend_covs)[,1]
trend_covs$lat <- coordinates(trend_covs)[,2]
trend_covs <- spTransform(trend_covs, CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
trend_covs <- trend_covs[,9:16]
results1 <- merge(res_sf, trend_covs@data, by.x="grid_id", by.y="amro_id", all=T)
results1[1310, 27:29] <- results1[1311, 27:29]
plot(results1["tempMean"])
names(results1)
cov_dat1 <- dplyr::select(results1, grid_id, na_id, x=centroid_x, y=centroid_y,
                   lon, lat, bcr_num=bcr, bcr_name, state,
                   country=country.x, alph, alph_ll, alph_ul, alph_iw,
                   eps, eps_ll, eps_ul, eps_iw, eps_sig, tau, tau_ll,
                   tau_ul, tau_iw, tau_sig, tempMean, popMean, geometry) %>%
  mutate(nLat=lat,
         tempMeanC=tempMean,
         popMeanP25=((popMean+1)/1000)^0.25,
         nLat_s=as.numeric(scale(lat)),
         tempMeanC_s=as.numeric(scale(tempMeanC)),
         popMeanP25_s=as.numeric(scale(popMeanP25)),
         car_int=grid_id) %>% arrange(grid_id)
# ------------------------------------------------------------------------------



# explore trend relationships --------------------------------------------------
# form 2
form2 <- tau ~ 1 + nLat_s + tempMeanC_s + popMeanP25_s +
  f(car_int, model="besag", graph=g1, scale.model=T, constr=T,
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01))))

# get results 2
out2<- inla(form2, family="gaussian", data=cov_dat1,
            control.compute=list(waic=T, config=T, cpo=T),
            control.inla=list(strategy="adaptive",
                              int.strategy="auto"),
            num.threads=3,
            verbose=F)
summary(out2)

# form 3
form3 <- tau ~ 1 + nLat + tempMeanC + popMeanP25 +
  f(car_int, model="besag", graph=g1, scale.model=T, constr=T,
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01))))

# get results 3
out3<- inla(form3, family="gaussian", data=cov_dat1,
            control.compute=list(waic=T, config=T, cpo=T),
            control.inla=list(strategy="adaptive",
                              int.strategy="auto"),
            num.threads=3,
            verbose=F)
summary(out3)
out3$waic$waic
-2*sum(out3$cpo$cpo, na.rm=T)
bri.hyperpar.summary(out3)
plot(out3)
hist(out3$cpo$pit)

# get random effects
cov_dat1$car_ints <- out3$summary.random$car_int$mean[1:n_cells]
cov_dat1$resid <- out3$summary.fitted.values$mean - cov_dat1$tau
cov_dat1$fit <- out3$summary.fitted.values$mean

# plot effects
plot(cov_dat1["car_ints"])
plot(cov_dat1["resid"])
plot(cov_dat1$fit, cov_dat1$tau)

# model fit
cov_dat1$fix_fit <- out3$summary.fixed$mean[1] +
  cov_dat1$nLat*out3$summary.fixed$mean[2] +
  cov_dat1$tempMeanC*out3$summary.fixed$mean[3] +
  cov_dat1$popMeanP25*out3$summary.fixed$mean[4]

ggplot(data=cov_dat1, aes(x=fit, y=tau)) +
  geom_point() + geom_abline(slope=1, intercept=0)
ggplot(data=cov_dat1, aes(x=fix_fit, y=tau)) +
  geom_point() + geom_abline(slope=1, intercept=0)

# r-squareds
cor(cov_dat1$tau, cov_dat1$fit, method="spearman", use="pairwise.complete.obs")
cor(cov_dat1$tau, cov_dat1$fix_fit, method="spearman", use="pairwise.complete.obs")

# posterior sampling function
fixed.posterior.sample <- function(result=out,
                                   full.ss=1000,
                                   part.ss=100){
  # set up
  out <- c()
  n.loops <- full.ss / part.ss
  all.par.names <- row.names(as.data.frame(
    inla.posterior.sample(1, result)[[1]]$latent))
  fix.names <- as.character(result$names.fixed)
  # get fix indices
  fix.idx <- c()
  for(i in fix.names){
    fix.idx <- c(fix.idx, grep(i, all.par.names))
  }
  # loop through samples
  for(i in 1:n.loops){
    post.1 <- inla.posterior.sample(part.ss, result)
    post.2 <- as.matrix(sapply(post.1, function(x) x$latent))
    fix.block <- t(post.2[min(fix.idx):max(fix.idx), ])
    colnames(fix.block) <- fix.names
    out <- rbind(out, fix.block)
    rm(post.1)
    rm(post.2)
    gc()
    print(paste("Subsample", i, "of", n.loops))
  }
  out <- as.data.frame(out)
  return(out)
}

# get fixed posteriors
posterior_ss <- 500
fixed_samp <- fixed.posterior.sample(result=out3,
                                     full.ss=posterior_ss,
                                     part.ss=100)
# get summaries of samples
apply(fixed_samp, 2, quantile, prob=c(0.025, 0.5, 0.975))

# make lat effects
fixed_samp <- as.matrix(fixed_samp)
newdat <- expand.grid(nLat=seq(min(cov_dat1$nLat, na.rm=T),
                                   max(cov_dat1$nLat, na.rm=T), length=50),
                      tempMeanC=mean(cov_dat1$tempMeanC, na.rm=T),
                      popMeanP25=mean(cov_dat1$popMeanP25, na.rm=T))
xmat <- model.matrix(~ nLat + tempMeanC + popMeanP25, data=newdat)
fitmat <- matrix(ncol=posterior_ss, nrow=nrow(newdat))
for(i in 1:posterior_ss) fitmat[,i] <- xmat %*% fixed_samp[i,]
newdat$med <- apply(fitmat, 1, quantile, prob=0.5)
newdat$lcl <- apply(fitmat, 1, quantile, prob=0.025)
newdat$ucl <- apply(fitmat, 1, quantile, prob=0.975)
p_lat <- ggplot(data=newdat, aes(x=nLat, y=med, ymin=lcl, ymax=ucl)) +
  geom_ribbon(fill="gray70") + geom_line()+
  xlab("Latititude (degrees north)") + ylab("Tau") +
  scale_y_continuous(limits=c(-7, 15)); p_lat

# make mean temp effects
newdat <- expand.grid(nLat=mean(cov_dat1$nLat, na.rm=T),
                      tempMeanC=seq(min(cov_dat1$tempMeanC, na.rm=T),
                                    max(cov_dat1$tempMeanC, na.rm=T), length=50),
                      popMeanP25=mean(cov_dat1$popMeanP25, na.rm=T))
xmat <- model.matrix(~ nLat + tempMeanC + popMeanP25, data=newdat)
fitmat <- matrix(ncol=posterior_ss, nrow=nrow(newdat))
for(i in 1:posterior_ss) fitmat[,i] <- xmat %*% fixed_samp[i,]
newdat$med <- apply(fitmat, 1, quantile, prob=0.5)
newdat$lcl <- apply(fitmat, 1, quantile, prob=0.025)
newdat$ucl <- apply(fitmat, 1, quantile, prob=0.975)
p_temp <- ggplot(data=newdat, aes(x=tempMeanC, y=med, ymin=lcl, ymax=ucl)) +
  geom_ribbon(fill="gray70") + geom_line()+
  xlab("Mean minimum winter temperature (C)") + ylab("Tau") +
  scale_y_continuous(limits=c(-7, 15)); p_temp

# make mean pop effects
newdat <- expand.grid(nLat=mean(cov_dat1$nLat, na.rm=T),
                      tempMeanC=mean(cov_dat1$tempMeanC, na.rm=T),
                      popMeanP25=seq(min(cov_dat1$popMeanP25, na.rm=T),
                                     max(cov_dat1$popMeanP25, na.rm=T), length=50))
xmat <- model.matrix(~ nLat + tempMeanC + popMeanP25, data=newdat)
fitmat <- matrix(ncol=posterior_ss, nrow=nrow(newdat))
for(i in 1:posterior_ss) fitmat[,i] <- xmat %*% fixed_samp[i,]
newdat$med <- apply(fitmat, 1, quantile, prob=0.5)
newdat$lcl <- apply(fitmat, 1, quantile, prob=0.025)
newdat$ucl <- apply(fitmat, 1, quantile, prob=0.975)
p_pop <- ggplot(data=newdat, aes(x=popMeanP25, y=med, ymin=lcl, ymax=ucl)) +
  geom_ribbon(fill="gray70") + geom_line()+
  xlab("Mean population size (power 0.025)") + ylab("Tau") +
  scale_y_continuous(limits=c(-7, 15)); p_pop

# plot model effects
multiplot(p_lat, p_temp, p_pop, cols=3)
pdf("trend_correlates.pdf", height=3.25, width=10.5)
multiplot(p_lat, p_temp, p_pop, cols=3)
dev.off()
# ------------------------------------------------------------------------------



# end analysis -----------------------------------------------------------------



