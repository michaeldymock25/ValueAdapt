
## this script produces the figures for the RSV case study (RSV.Rmd) (note: the code in RSV.Rmd must be run first)

library(dplyr)

set.seed(12856421)

n <- seq(50, 10000, by = 50)
out_perfect <- vector(length = length(n))
out_sample <- vector(length = length(n))
for(i in 1:length(n)){
  print(n[i])
  out_perfect[i] <- enb_perfect(D = D, U = U, Theta = prior, t = t[1,], prop = prop, cost = delta + gamma*n[i])
  samp_fun <- function(args){
    c(II = rbinom(1, size = n[i], prob = args[["p_II"]]),
      MV = rbinom(1, size = n[i], prob = args[["p_MV"]]))
  }
  out_sample[i] <- enb_sample(D = D, U = U, Theta = prior, t = t[1,], prop = prop, cost = delta + gamma*n[i],
                              samp_fun = samp_fun, stat_fun = stat_fun, model = model)
}
dat_vis <- data.frame(n = rep(n, 2),
                      value = c(out_sample, out_perfect),
                      grp = rep(c("ENBS", "ENBP"), each = length(n)))

tiff("Examples/RSV/Figure 2.tiff", width = 6, height = 7, units = "in", res = 1000)

ggplot(dat_vis, aes(x = n, y = value, linetype = grp)) +
  geom_line() +
  scale_x_continuous("Proposed interim sample size per arm", n.breaks = 5) +
  scale_y_continuous("ENB", n.breaks = 5) +
  scale_linetype_manual("", values = c("dashed", "solid")) +
  theme_bw() +
  theme(legend.position = "bottom")

dev.off()

post_vis <- data.frame(Analysis = 0,
                       value = prior[,"p_II"] - prior[,"p_MV"])
post_vis_1 <- data.frame(Scenario = "Scenario 1",
                         Analysis = rep(1:4, each = nrow(post_1_long)/4/2),
                         value = post_1_long[post_1_long$Parameter == "p_MV",]$x -
                                 post_1_long[post_1_long$Parameter == "p_II",]$x)
post_vis_2 <- data.frame(Scenario = "Scenario 2",
                         Analysis = rep(1:4, each = nrow(post_2_long)/4/2),
                         value = post_2_long[post_2_long$Parameter == "p_MV",]$x -
                                 post_2_long[post_2_long$Parameter == "p_II",]$x)
post_vis <- bind_rows(post_vis |> mutate(Scenario = "Scenario 1"),
                      post_vis |> mutate(Scenario = "Scenario 2"),
                      post_vis_1,
                      post_vis_2)
post_vis$Analysis <- factor(post_vis$Analysis, labels = c("Prior", "Analysis 1", "Analysis 2", "Analysis 3", "Analysis 4"))
post_vis$int <- ifelse(post_vis$Scenario == "Scenario 1", 0.08, 0.02)

tiff("Examples/RSV/Figure 3.tiff", width = 6, height = 7, units = "in", res = 1000)

ggplot(post_vis, aes(x = value)) +
  facet_grid(Analysis ~ Scenario) +
  geom_density() +
  geom_vline(aes(xintercept = int), linetype = "dashed") +
  scale_x_continuous("Distribution of the absolute risk difference", limits = c(-0.2, 0.2)) +
  ylab("Density") +
  theme_bw() +
  theme(legend.position = "bottom")

dev.off()

net_benefit_vis <- data.frame(Analysis = "Prior", INMB = INMB_0, E = mean(INMB_0))
net_benefit_vis <- bind_rows(net_benefit_vis |> mutate(Scenario = "Scenario 1"),
                             net_benefit_vis |> mutate(Scenario = "Scenario 2"),
                             INMB_1 |> mutate(Scenario = "Scenario 1") |> rename(INMB = x),
                             INMB_2 |> mutate(Scenario = "Scenario 2") |> rename(INMB = x))
net_benefit_vis$Analysis <- factor(net_benefit_vis$Analysis, levels = c("Prior", "Analysis 1", "Analysis 2", "Analysis 3", "Analysis 4"))

tiff("Examples/RSV/Figure 4.tiff", width = 6, height = 7, units = "in", res = 1000)

ggplot(net_benefit_vis, aes(x = INMB)) +
  facet_grid(Analysis ~ Scenario) +
  geom_density() +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = E), linetype = "dashed", colour = "#CC6677") +
  scale_x_continuous("INMB", limits = c(-2500, 2500)) +
  ylab("Density") +
  theme_bw() +
  theme(legend.position = "bottom")

dev.off()

tiff("Examples/RSV/Figure 5.tiff", width = 6, height = 7, units = "in", res = 1000)

ggplot(sims_prop, aes(x = Analysis, y = Prop, group = Scenario, colour = Scenario)) +
  geom_line() +
  scale_y_continuous("Proportion of trials stopped", limits = c(0,1)) +
  scale_colour_manual("", values = colours) +
  theme_bw() +
  theme(legend.position = "bottom")

dev.off()

tiff("Examples/RSV/Figure 6.tiff", width = 6, height = 7, units = "in", res = 1000)

ggplot(sims_summ, aes(x = Analysis, y = mean, ymin = lower, ymax = upper, colour = Measure, group = Measure)) +
  facet_wrap(~Scenario) +
  geom_pointrange(position = position_dodge(width = 0.25)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual("", values = colours) +
  scale_y_continuous("Value (1 million AUD)", n.breaks = 6) +
  theme_bw() +
  theme(legend.position = "bottom")

dev.off()