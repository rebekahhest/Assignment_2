##Assignment 1

#packages -----
library(readr)
library(dplyr)
library(tidyr)
library(vegan)
library(stringr)
library(mgcv)
library(ggplot2)

####retrieve public data via BOLD API ------
#retrieved on October 9, 2025
#parasitoid dataset
  #dfIchneumoninae <- read_tsv("https://www.boldsystems.org/index.php/API_Public/specimen?taxon=Ichneumoninae&format=tsv")
  #write_tsv(dfIchneumoninae, "../data/df_Ichneumoninae.tsv")

df_Ichneumoninae <- read_tsv("../data/df_Ichneumoninae.tsv")

#host dataset
  #dfHeliconiinae_Danainae <- read_tsv("https://www.boldsystems.org/index.php/API_Public/specimen?taxon=Heliconiinae|Danainae&format=tsv")
  #write_tsv(dfHeliconiinae_Danainae, "../data/df_Heliconiinae_Danainae.tsv")

df_Heliconiinae_Danainae <- read_tsv("../data/df_Heliconiinae_Danainae.tsv")

####filter data -----
parasitoid <- df_Ichneumoninae %>%
  select("processid", "sampleid", "bin_uri", "genus_name", "lat", "lon") %>%
  filter(!is.na(bin_uri) & !is.na(lat) & !is.na(lon))
write_tsv(parasitoid, "../data/parasitoid.tsv")

host <- df_Heliconiinae_Danainae %>%
  select("processid", "sampleid", "bin_uri", "genus_name", "lat", "lon") %>%
  filter(!is.na(bin_uri) & !is.na(lat) & !is.na(lon))
write_tsv(host, "../data/host.tsv")

rm(df_Ichneumoninae, df_Heliconiinae_Danainae)

#latitude histogram to check sampling bias between datasets
hist(parasitoid$lat, main = "Sampling Distribution for Parasitoid Dataset", xlab = "Latitude (°)")
hist(host$lat, main = "Sampling Distribution for Host Dataset", xlab = "Latitude (°)")
#similar distribution so sampling bias is less of a concern

#####create 10 degree latitude bins where bounds are nearest multiple of 10 ------

##find lower and upper bin bounds for all data
#lower bound of lowest lat bin
  #1) include lowest lat value: a = min(lat)
  #2) transform coordinate system into multiple of 10: b = a / 10
  #3) use nearest lower integer (captures decimals): c = floor(b)
  #4) convert back to coordinate scale: c * 10
parasitoid_min_lat <- floor(min(parasitoid$lat) / 10) * 10
host_min_lat <- floor(min(host$lat) / 10) * 10
#upper bound of highest lat bin
  #1) include highest lat value: a = max(lat)
  #2) transform coordinate system into multiple of 10: b = a / 10
  #3) use nearest higher integer (captures decimals): c = ceiling(b)
  #4) convert back to coordinate scale: c * 10
parasitoid_max_lat <- ceiling(max(parasitoid$lat) / 10) * 10
host_max_lat <- ceiling(max(host$lat) / 10) * 10

#overall range of data spans from [-50, 90] degrees

#creates bins using the lower and upper bound for all data where each interval increment is 10
parasitoid_intervals <- seq(parasitoid_min_lat, parasitoid_max_lat, by = 10)
host_intervals <- seq(host_min_lat, host_max_lat, by = 10)

rm(parasitoid_min_lat, host_min_lat, parasitoid_max_lat, host_max_lat)

#organize processids into latitude bins
parasitoid_lat_bins <- parasitoid %>%
  mutate(
    lat_bin = cut(
      lat,
      breaks = parasitoid_intervals,
      include.lowest = TRUE,
      right = FALSE
    )
  )

host_lat_bins <- host %>%
  mutate(
    lat_bin = cut(
      lat,
      breaks = host_intervals,
      include.lowest = TRUE,
      right = FALSE
    )
  )

rm(parasitoid_intervals, host_intervals)

####calculate Shannon diversity ------
#structure data frame 
parasitoid_matrix <- parasitoid_lat_bins %>%
  group_by(lat_bin, bin_uri) %>%
  summarise(count = n()) %>% 
  pivot_wider(names_from = bin_uri, values_from = count, values_fill = 0)

host_matrix <- host_lat_bins %>%
  group_by(lat_bin, bin_uri) %>%
  summarise(count = n()) %>% 
  pivot_wider(names_from = bin_uri, values_from = count, values_fill = 0)

rm(parasitoid_lat_bins, host_lat_bins)

#calculate Shannon Index
parasitoid_div <- diversity(x = parasitoid_matrix[,-1], index = "shannon")
host_div <- diversity(x = host_matrix[,-1], index = "shannon")

#create data frame to only include latitude bins and shannon index
parasitoid_shannon <- data.frame(
  lat_bin = parasitoid_matrix$lat_bin,
  shannon = parasitoid_div
)

host_shannon <- data.frame(
  lat_bin = host_matrix$lat_bin,
  shannon = host_div
)

rm(parasitoid_matrix, parasitoid_div, host_matrix, host_div)

#transform lat_bins into one numeric x value for regression
#midpoint used as representative x value for lat_bin single value
parasitoid_variables <- parasitoid_shannon %>% 
  rowwise() %>%
  mutate(min = as.numeric(unlist(str_extract_all(lat_bin, "-?\\d+"))[1]),
         max = as.numeric(unlist(str_extract_all(lat_bin, "-?\\d+"))[2]),
         lat_midpoint = mean(c(min,max))) %>%
  select(lat_midpoint, shannon)

host_variables <- host_shannon %>% 
  rowwise() %>%
  mutate(min = as.numeric(unlist(str_extract_all(lat_bin, "-?\\d+"))[1]),
         max = as.numeric(unlist(str_extract_all(lat_bin, "-?\\d+"))[2]),
         lat_midpoint = mean(c(min,max))) %>%
  select(lat_midpoint, shannon)

####Generalized Additive Models (GAMs) ------
#parasitoid diversity as predicted by latitude
parasitoid_lat_gam <- gam(shannon ~ s(lat_midpoint), data = parasitoid_variables, family = gaussian())
summary(parasitoid_lat_gam)
#p-value 0.0651 > 0.05; not significant

#host diversity as predicted by latitude
host_lat_gam <- gam(shannon ~ s(lat_midpoint), data = host_variables, family = gaussian())
summary(host_lat_gam)
#p-value 0.00153 < 0.05; significant

#merge datasets grouped by common column latitude
merged_lat_shannon <- merge(parasitoid_variables, host_variables, by = "lat_midpoint") %>%
  rename(parasitoid = shannon.x) %>%
  rename(host = shannon.y) %>%
  pivot_longer(cols = c(parasitoid, host),
               names_to = "group", values_to = "shannon")

write_tsv(merged_lat_shannon, "../data/merged_lat_shannon.tsv")

#ensure groups are ordered by parasitoid then host, not alphabetically (default)
merged_lat_shannon$group <- factor(merged_lat_shannon$group,
                                   levels = c("parasitoid", "host"))

#plot both datasets where shannon index is predicted by latitude
ggplot(merged_lat_shannon, aes(x = lat_midpoint, y = shannon, colour = group, shape = group)) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey") +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = c("parasitoid" = "red", "host" = "blue"),
                      labels = c("Parasitoid", "Host"),
                      name = "Dataset") +
  scale_shape_manual(values = c("parasitoid" = 17, "host" = 16),
                     labels = c("Parasitoid", "Host"),
                      name = "Dataset") +
  labs(title = "Diversity Across Latitudes",
       x = "Latitude (°)",
       y = "Shannon Index (H')",
       colour = "Dataset") + 
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("../figs/fig1_GAM1_GAM2.png", plot = last_plot() , width = 12, height = 6, dpi = 100)
  
#merge datasets for generalizd additive model without pivot_longer
merged_gam <- merge(parasitoid_variables, host_variables, by = "lat_midpoint") %>%
  rename(parasitoid_shannon = shannon.x) %>%
  rename(host_shannon = shannon.y)

write_tsv(merged_gam, "../data/merged_gam.tsv")

####Generalized Additive Model
#parasitoid diversity as predicted by latitude and host diversity
gam_model <- gam(parasitoid_shannon ~ s(lat_midpoint) + s(host_shannon), data = merged_gam)
summary(gam_model)
#s(lat_midpoint) p-value 0.0474 < 0.05; significant
#s(host_shannon) p-value 0.0231 < 0.5; significant

# Plot the partial effect of each independent variable

#latitude
plot(gam_model, select = 1, 
     xlab = "Latitude (°)", 
     ylab = "Partial Effect on Parasitoid Diversity")
png("../figs/fig2_GAM3_latitude.png", width = 1200, height = 900, res = 150)
plot(gam_model, select = 1, 
     xlab = "Latitude (°)", 
     ylab = "Partial Effect on Parasitoid Diversity")
dev.off()


#host_shannon
plot(gam_model, select = 2,
     xlab = "Host Shannon Index (H')", 
     ylab = "Partial Effect on Parasitoid Diversity")
png("../figs/fig3_GAM3_host_div.png", width = 1200, height = 900, res = 150)
plot(gam_model, select = 2,
     xlab = "Host Shannon Index (H')", 
     ylab = "Partial Effect on Parasitoid Diversity")
dev.off()