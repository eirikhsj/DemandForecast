

#PCA-codes


#Alternative 1: Do SVD straing on x

prep = prep_demand_temp_data(include_na_volume = FALSE)
#
X_mat = prep$X_mat
X_train = X_mat[I_train,]

X_train = X_mat
##--- Center the data ---
mu = colMeans(X_train)

X_center = matrix(NA, dim(X_train)[1], dim(X_train)[2])
for(j in 1:dim(X_train)[2]){
    X_center[, j] = X_train[, j] - mu[j]
}
##---- Perform the SVD ---
a = Sys.time();

l_svd = svd(X_center)# 59.53304 sec
# l_svd = svd(X_center)
U = l_svd$v[, 1]
print(Sys.time()-a)

F_1  = X_center %*% U
plot(F_1, type = 'l')

cumsum(l_svd$d^2/sum(l_svd$d^2))[1:10]
(l_svd$d^2/sum(l_svd$d^2))[1:10]


cumsum(l_svd$d^2/sum(l_svd$d^2))[1:10]
(l_svd$d^2/sum(l_svd$d^2))[1:10]


#Alternative 2: Perform SVD on Sigma
a = Sys.time();

Sigma = t(X_center) %*% X_center
l_svd_sigma = svd(Sigma)
print(Sys.time()-a) #16.64628
U_sig = l_svd_sigma$u[, 1]

F_1_sig  = X_center %*% U_sig

cumsum(l_svd_sigma$d/sum(l_svd_sigma$d))[1:10]
(l_svd_sigma$d/sum(l_svd_sigma$d))[1:10]

singular_values <- l_svd_sigma$d[1:15]/(length(X_center[,1])-1)
total_variance <- sum(l_svd_sigma$d/(length(X_center[,1])-1))
total_explained <- cumsum(singular_values / total_variance)
variance_explained <- singular_values / total_variance
# Create a data frame for the scree plot

scree_data <- data.frame(Component = 1:length(singular_values),
                         Variance = singular_values,
                         Total = total_explained,
                         Explained = variance_explained)

scree_data$Explained_Scaled <- scree_data$Explained * total_variance
scree_data$Total_Scaled <- scree_data$Total * total_variance
scree_data$Total[1] = -100000
scree_data$Total_Scaled[1] = -100000
# Create the scree plot using ggplot2
pdf("/Users/Eirik/Desktop/Master2023/ExportedPlots/Scree_plot.pdf", height = 7, width = 9, pointsize = 12)

ggplot(scree_data) +

    geom_bar(aes(x = Component, y = Total_Scaled), stat = "identity", fill = "grey", alpha = 0.5) +
    geom_bar(aes(x = Component, y = Explained_Scaled), stat = "identity", fill = "red", alpha = 0.5) +
    geom_line(aes(x = Component, y = Variance), color = "blue") +
    geom_point(aes(x = Component, y = Variance), color = "blue") +
    geom_text(aes(x = Component, y = Explained_Scaled, label = scales::percent(Explained, accuracy = 0.1)),
              vjust = -0.8, size = 3) +
    geom_text(aes(x = Component, y = Total_Scaled, label = scales::percent(Total, accuracy = 0.1)),
              vjust = -0.8, size = 3) +
    labs(x = "Principal Components", y = "Variance", title = "Scree Plot of In-Sample ERA5 Temperature Principal Components") +
    scale_x_continuous(breaks = 1:max(scree_data$Component),
                       labels = function(x) ifelse(x %in% c(1, 5, 10, 15), x, "")) +
    scale_y_continuous(limits = c(0, total_variance+1000),sec.axis = sec_axis(trans = ~ .,
                                                                              name = "Explained Variance (%)",
                                                                              breaks = total_variance*c(0, 0.2, 0.4, 0.6, 0.8, 1),labels = c(0, 20, 40, 60, 80, 100))) +
    theme_bw()

dev.off()

#Alternative 3: Perform SVD on Sigma
Sigma2 = (t(X_center) %*% X_center) /(length(X_center[,1])-1)
l_svd_sigma2 = svd(Sigma2)
U_sig2 = l_svd_sigma2$u[, 1]

cumsum(l_svd_sigma2$d/sum(l_svd_sigma2$d))[1:10]
(l_svd_sigma2$d/sum(l_svd_sigma2$d))[1:10]

F_1_sig2  = X_center %*% U_sig2

F_1_sig2

#Alt 4 incorporate n-1
Sigma = t(X_center) %*% X_center
l_svd_sigma = svd(Sigma)
print(Sys.time()-a) #16.64628
U_sig = l_svd_sigma$u[, 1]

F_1_sig  = X_center %*% U_sig

cumsum(l_svd_sigma$d/sum(l_svd_sigma$d))[1:10]
(l_svd_sigma$d/sum(l_svd_sigma$d))[1:10]


hist(l_svd$v[, 1])


### GAM plotting
all_temp_data = get_pca(X_mat, seq(1:nrow(X_mat)), p_comps = 460)
plotgam = gam(volume~ s(PC1), data =all_temp_data$dt_train)
plotgam2 = gam(volume~ s(PC2), data =all_temp_data$dt_train)

plotgam = getViz(plotgam)
o <- plot( sm(plotgam, 1) )
o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

qq(plotgam, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2),
   ylim = c(-25000, 25000),
   xlim = c(-25000, 25000))


check(plotgam,
      a.qq = list(method = "tnorm",
                  a.cipoly = list(fill = "light blue")),
      a.respoi = list(size = 0.2, col = 'blue'),
      a.hist = list(bins = 30, fill = "blue", color = "grey", title = "Histogram of Residu"))


plotgam2 = getViz(plotgam2)
o <- plot( sm(plotgam2, 1) )
o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

load("~/Desktop/Master2023/Output_data/New_Best_NWP_model_output/w_n.Rda")
weight_plots = data.table()
for(i in i:24){
    setorder(weight_plot, by = 'V1')
    namez = paste0('day', i)
    weight_plots[,(namez) := weight_plot$weight ]
}

library(ggdist)
# Parameters of the beta distribution
alpha <- 2
beta <- 5

# Generate data for the PDF
x <- seq(0, 1, length.out = 100)
pdf <- dbeta(x, alpha, beta)

# Create a data frame for plotting
df <- data.frame(x = x, pdf = pdf)

# Create the plot
ggplot(df, aes(x = x, y = pdf)) +
    geom_line() +
    geom_area(data = subset(df, x >= 0.3 & x <= 0.5),
              aes(fill = "Green"), alpha = 0.5) +
    scale_fill_manual(values = c("Green" = "green"),
                      guide = guide_none()) +
    labs(x = "x", y = "Density", title = "PDF of Beta Distribution") +
    theme_minimal()

acf(resid(plotgam$lme))
plot(gam.check(plotgam)$residuals, gam.check(plotgam)$fitted.values, xlab = "Fitted values", ylab = "Residuals")

plotgam = gam(volume~ s(PC1, k = 20), data =all_temp_data$dt_train) #Check with more
gam.check(plotgam)
appraise(plotgam)


df <- data.frame(
    month = rep(month.abb, each = 100),
    value = rnorm(1200)
)

# create the plot
ggplot(df, aes(x = month, y = value, fill = month)) +
    geom_violin() +
    facet_wrap(~ month, ncol = 3, scales = "free_x") +
    theme_minimal()


library("ggplot2")
library("dplyr")
sm <- smooth_estimates(plotgam)
sm %>%
    add_confint() %>%
    ggplot(aes(y = est, x = x2)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
                alpha = 0.2, fill = "forestgreen") +
    geom_line(colour = "forestgreen", size = 1.5) +
    labs(y = "Partial effect",
         title = expression("Partial effect of" ~ f(x[1])),
         x = expression(x[1]))

mat <- matrix(U2, nrow = 21, ncol = 22, byrow = TRUE)
ggplot(melt(mat), aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue")

plot(matrix(l_svd$u, nrow = 462, ncol = 1, byrow = TRUE), type = 'l')
lines(matrix(l_svd$v, nrow = 462, ncol = 1, byrow = TRUE), col = 'red')
plot(matrix(l_svd2$u, ncol = 1, byrow = FALSE))


qqline(date_demand$volume, rt(., 10), df = 10)
summary(date_demand$volume)
shap = shapiro.test(date_demand$volume[1:5000])
ks.test(date_demand$volume, "pnorm")

SSS = cor(X_mat)
heatmap(SSS, col = heat.colors(256), Rowv = NA, Colv = NA)

df <- reshape2::melt(SSS)

# Create a heatmap using ggplot2
ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = heat.colors(256))

ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "yellow")+
    labels(ylab = "Grid Point")


#
acf_values <- apply(mat, 2, acf)

# Access the autocorrelation values for the first column
acf_values[[2]]$acf


#Histogram demand
ggplot(data.frame(x = vector), aes(x)) +
    geom_density(date_demand[month(date) == 1,volume])
abline(v = mean(date_demand$volume), col = 'red')
abline(v = mode(date_demand$volume), col = 'red')

df <- data.frame(
    value = c(vector1, vector2, vector3),
    group = rep(c("Vector 1", "Vector 2", "Vector 3"), each = 1000)
)

# Create density plot with overlapping curves
ggplot(df, aes(x = value, color = group)) +
    geom_density() +
    scale_color_manual(values = c("red", "blue", "green"))

ggplot(df, aes(x = value, color = group, fill = group)) +
    geom_density(alpha = 0.5) +
    scale_color_manual(values = c("red", "blue", "green")) +
    scale_fill_manual(values = c("red", "blue", "green"))


data1 <- date_demand[month(date) %in% c(12,1,2), volume]
data2 <- date_demand[month(date) %in% c(3,4,5), volume]
data3 <- date_demand[month(date) %in% c(6,7,8), volume]
data4 <- date_demand[month(date) %in% c(9,10,11), volume]
data5 <- date_demand[, volume]

# Determine the maximum length among the vectors
max_length <- length(data5)

# Pad the shorter vectors with NA values
data1 <- c(data1, rep(NA, max_length - length(data1)))
data2 <- c(data2, rep(NA, max_length - length(data2)))
data3 <- c(data3, rep(NA, max_length - length(data3)))
data3 <- c(data3, rep(NA, max_length - length(data4)))
data3 <- c(data3, rep(NA, max_length - length(data5)))
# Combine data into a single data frame
combined_data <- data.frame(Value = c(data1, data2, data3, data4, data5),
                            Period = rep(c("Winter", "Spring", "Summer","Fall","Overall"),
                                        each = max_length))

# alpha_values <- c("Group 1" = 0.2, "Group 2" = 0.3, "Group 3" = 0.4, "Group 4" = 0.5, "Group 5" = 1)
colors <- c("Winter" = "blue", "Spring" = "green", "Summer" = "red","Fall" = "brown","Overall" = "grey")
# Plotting
pdf("/Users/Eirik/Desktop/Master2023/ExportedPlots/ElectricityDensity.pdf", height = 7, width = 7, pointsize = 24)

ggplot(combined_data, aes(x = Value,y = ..count.., fill = Period)) +
    geom_density(alpha = c(0.3)) +
    labs(title = "Electrity Demand Density Plot",
         x = "Electricity Demand (MWh)", y = "Density")+
    scale_fill_manual(values = colors)+
    theme(panel.background = element_blank())+

    theme(legend.position = c(0.8,0.8))

    theme_minimal()
dev.off()


pdf("/Users/Eirik/Desktop/Master2023/ExportedPlots/ElectricityDemand.pdf", height = 7, width = 7,
    pointsize = 24)
# All points plot wiht line
    ggplot(data.frame(x = 1:dim(date_demand)[1], y = date_demand$volume), aes(x = x, y = y)) +
        geom_point(size = 0.3, color = "lightblue", alpha = 0.2,show.legend = TRUE, legend_key = TRUE) +
        geom_line(data = data.frame(y = c(predict(reg_dem),rep(NA,744)),x = 1:dim(date_demand)[1]), aes(x,y), col = 'red')+
        labs(title = "Norwegian Electrity Demand (Jan 2013 - Jan 2023)",
             x = "Year", y = "Electricity Demand (MWh)")+
        scale_x_continuous(breaks = seq(1,88392,8760), labels = seq(2013,2023, by = 1))+
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank())

dev.off()


pdf("/Users/Eirik/Desktop/Master2023/ExportedPlots/ACF_NORD.pdf", height = 7, width = 9, pointsize = 12)
ggacf_plot_ENER = ggAcf(diff(ts(date_demand$volume,frequency = 1)), lag.max = 500)+
    labs(title = "Autocorrelation Function Energy Demand",
         x = "Lags (1-hour intervals)",
         y = "ACF")
#theme_bw() +
#scale_x_continuous(breaks=seq(0, 60, by = 10)*4,labels=c(0,5,10,15, 20,25,30))+

# theme(plot.title = element_text(hjust = 0.5),
#       axis.line = element_line(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank())
ggacf_plot_ENER
dev.off()


winter = date_demand[month(date)==2,mean(volume), by = hour]
summer = date_demand[month(date)==6,mean(volume), by = hour]
fall = date_demand[month(date)==10,mean(volume), by = hour]
spring = date_demand[month(date)==4,mean(volume), by = hour]
max_length_hour = 24
combined_data_hour <- data.frame(x = 1:24, winter = winter$V1, spring = spring$V1,
                                 summer = summer$V1, fall= fall$V1)

ggplot(combined_data_hour, aes(x = x, y = winter)) +
    geom_line(color = 'blue')+
    geom_line(aes(y = spring), color = "green") +
    geom_line(aes(y = summer), color = "red") +
    geom_line(aes(y = fall), color = "brown") +
    scale_x_continuous(breaks = c(4,8,12,16,20,24), labels = c(4,8,12,16,20,24))+

    labs(title = "Mean Electrity Demand by Hour (Jan 2013 - Jan 2023)",
         x = "Hour", y = "Electricity Demand (MWh)")+
    theme(legend.position = c(0.8,0.8))+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())


winter = date_demand[month(date)==2,mean(volume), by = wday(date)]
summer = date_demand[month(date)==6,mean(volume), by = wday(date)]

fall = date_demand[month(date)==10,mean(volume), by = wday(date)]
spring = date_demand[month(date)==4,mean(volume), by = wday(date)]

setnames(winter, old = "V1", new = "winter")
setnames(spring, old = "V1", new = "spring")

setnames(fall, old = "V1", new = "fall")

setnames(summer, old = "V1", new = "summer")

merged_dt <- merge(winter, summer, by = "wday")
merged_dt <- merge(merged_dt, fall, by = "wday")
merged_dt <- merge(merged_dt, spring, by = "wday")

merged_dt$x = merged_dt$wday-1
merged_dt[x==0,x:= 7]
combined_data_week<- data.frame(merged_dt[2:6])

ggplot(combined_data_hour, aes(x = x, y = winter)) +
    geom_line(color = 'blue')+
    geom_line(aes(x = x, y = spring), color = "green") +
    geom_line(aes(x = x, y = summer), color = "red") +
    geom_line(aes(x = x, y = fall), color = "brown") +
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c(1,2,3,4,5,6,7))+

    labs(title = "Mean Electrity Demand by Hour (Jan 2013 - Jan 2023)",
         x = "Weekday", y = "Electricity Demand (MWh)")+
    theme(legend.position = c(0.8,0.8))+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())








# Extract year and month from date column
ERA_NWP_RQE1$year <- format(ERA_NWP_RQE1$date, "%Y")
ERA_NWP_RQE1$month <- format(ERA_NWP_RQE1$date, "%m")

# Create a factor variable for month with correct order
ERA_NWP_RQE1$month <- factor(ERA_NWP_RQE1$month, levels = sprintf("%02d", 1:12))

# Compute mean value for each month
mean_data <- aggregate(mean_loss_1_19 ~ year + month, ERA_NWP_RQE1, mean)

# Create bar plot
ggplot(mean_data, aes(x = month, y = mean_loss_1_19, fill = year)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_discrete(name = "Year") +
    labs(x = "Month", y = "Mean Value")
