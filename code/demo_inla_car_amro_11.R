# cbc with car model -----------------------------------------------------------
# this is code for running an ICAR SVC model on CBC data.
# ------------------------------------------------------------------------------



# set up -----------------------------------------------------------------------
# install the required packages
library(tidyverse)
library(RColorBrewer)
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
cells_with_counts <- unique(modeling_data$grid_id[which(
  !is.na(modeling_data$hours))])
alph <- exp(out1$summary.random$alpha_i$`0.5quant`[cells_with_counts])
alph_ll <- exp(out1$summary.random$alpha_i$`0.025quant`[cells_with_counts])
alph_ul <- exp(out1$summary.random$alpha_i$`0.975quant`[cells_with_counts])
alph_iw <- alph_ul - alph_ll
hist(alph); summary(alph)
eps <- out1$summary.random$eps_i$`0.5quant`[cells_with_counts]
eps_ll <- out1$summary.random$eps_i$`0.025quant`[cells_with_counts]
eps_ul <- out1$summary.random$eps_i$`0.975quant`[cells_with_counts]
eps_iw <- eps_ul - eps_ll
hist(eps); summary(eps)
hist(eps_ll); summary(eps_ll)
hist(eps_ul); summary(eps_ul)
tau <- (exp(out1$summary.random$tau_i$`0.5quant`[cells_with_counts])
        - 1) * 100
tau_ll <- (exp(out1$summary.random$tau_i$`0.025quant`[cells_with_counts])
           - 1) * 100
tau_ul <- (exp(out1$summary.random$tau_i$`0.975quant`[cells_with_counts])
           - 1) * 100
tau_iw <- tau_ul - tau_ll
hist(tau); summary(tau)
hist(tau_ll); summary(tau_ll)
hist(tau_ul); summary(tau_ul)
hist(tau_iw); summary(tau_iw)

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
#     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
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
# out2$summary.random$grid_yr[grep("600-", out2$summary.random$grid_yr$ID,
#                                  fixed=T), c(1,5,4,6)]
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
results_cells <- merge(cbc_amro_grid, post_sum)
res_sf <- as(results_cells, "sf")

# filter out empty cells
cells_with_counts
#hits1 <- which(!is.na(over(cbc_amro_grid, sd2)["grid_id"])==T)
cbc_amro_grid <- cbc_amro_grid[cells_with_counts, ]
plot(cbc_amro_grid)
plot(cbc_amro_circles, add=T, pch=16, cex=0.5)
res_sf <- res_sf[cells_with_counts, ]
summary(post_sum)

# map tau
tau_p1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau), col="gray40", size=0.3) +
  scale_fill_gradient2("Tau\n(% per year)", low = ("royalblue4"),
                       mid = "white",
                       high = ("red4"), midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau_p2 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_iw), col="gray40", size=0.3) +
  scale_fill_gradient2("Tau\ncredible\ninterval\nwidth\n(% per year)",
                       low = ("purple4"), mid = "white",
                       high = ("green4"), midpoint = 6, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
tau_p3 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=tau_sig), col="gray40", size=0.3) +
  scale_fill_gradient2("Significant tau\n(% per year)",
                       low = muted("royalblue4"), mid = "gray95",
                       high = muted("red4"), midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar") +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
multiplot(tau_p1)
# multiplot(tau_p1, tau_p3, tau_p2, cols=2)

# map epsilon
eps_p1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=eps), col="gray40", size=0.3) +
  scale_fill_gradient2("Epsilon", low = muted("purple4"), mid = "white",
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
multiplot(eps_p1)
# multiplot(eps_p1, eps_p3, eps_p2, cols=2)

# map alpha
alph_p1 <- ggplot() +
  geom_sf(data=bcr_sf, fill="gray40", col="gray40") +
  geom_sf(data=res_sf, aes(fill=alph), col="gray40", size=0.3) +
  scale_fill_gradient2("Alpha", low = "tan4", mid = "white",
                       high = "green4", midpoint = 0, space = "Lab",
                       na.value = "grey40", guide = "colourbar", trans="log",
                       breaks=c(0.01, 0.1, 1, 10)) +
  theme_map() + theme(panel.grid.major=element_line(colour="transparent"))
multiplot(alph_p1)

# print cell maps
pdf("tau_alpha.pdf", height=7.5, width=10.5)
multiplot(tau_p1, tau_p3, tau_p2, alph_p1, cols=2)
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
  select(1,5,6,7,11) %>%
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



# end analysis -----------------------------------------------------------------



