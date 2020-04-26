######################
## Import Libraries ##
######################

library(readxl)
library(tidyverse)
library(lubridate)
library(nlme)
library(grid)
library(gridExtra)
library(ggpmisc)

####################
## Importing Data ##
####################

#importing original data; set to correct working directory
chicken_orig <- read_excel("./Air movement temperature data.xlsx", 
                           sheet = "Temperature data")  

#############################
## Variable Transformation ##
#############################

#this is the 13,360-observation data-set for flank-rectal investigation
chicken_trans <-
  chicken_orig %>% 
  
  #renaming variables
  select(day = `Day of age`, time_of_day = `date-time`, 
         body_temp = Temperature, set = Set, 
         loc = `location R = rectal, W = wing`,
         sensor =`Sensor number`, room_temp = `Room temperature (F)`, 
         trial = Trial, rep = Rep, parent_flock = `Parent flock age`,
         group = Group, period = `Period ( P=pen; T=treatment period)`,
         air_velocity = `Air velocity (ft/min)`
  ) %>%
  
  #recoding, creating chick id variable
  mutate(time_of_day = if_else(hour(time_of_day) < 12, 'AM', 'PM'),
         id = paste(paste(sensor, rep, sep = ""), 
                    room_temp, sep = "")) %>% 
  
  #creating factor variables
  mutate_at(c("day", "time_of_day", "set", "loc", 
              "sensor", "room_temp", "trial", "rep", 
              "parent_flock", "group", "period", 
              "air_velocity"), as_factor) %>%
  
  #placing flank and rectal temperatures on same line
  pivot_wider(names_from = loc, values_from = body_temp, 
              names_prefix = "body_temp_") %>%
  
  #creating flank - rectal difference (adj is rounded to 0.1)
  mutate(body_temp_diff = body_temp_W - body_temp_R,
         body_temp_R_adj = round(body_temp_R, 1),
         body_temp_W_adj = round(body_temp_W, 1)) %>%  
  mutate(body_temp_diff_adj = round(body_temp_W_adj 
                                    - body_temp_R_adj, 1)) %>%
  
  #factor version of body temp diff
  mutate(body_temp_diff_fac = factor(body_temp_diff_adj))

#this is the 2,672-observation data for modeling cell body temperature
chicken_agg <- 
  chicken_trans %>%
  
  #dropping unneeded variables
  select(-c(body_temp_diff, body_temp_diff_adj, body_temp_diff_fac,
            body_temp_R_adj, body_temp_W_adj,
            trial, air_velocity, sensor)) %>%
  
  #aggregating data
  group_by(rep, day, time_of_day, room_temp, group, period, id) %>%
  summarize_at(c("body_temp_W", "body_temp_R"), mean) %>%
  pivot_longer(starts_with("body_temp"), 
               names_to = "loc", 
               values_to = "body_temp") %>%
  pivot_wider(names_from = period, 
              values_from = body_temp, 
              names_prefix = "body_temp_") %>%
  ungroup %>%
  
  #placing flank and rectal measurements on different lines
  mutate_at("loc", as_factor) %>%
  
  #recoding variables
  mutate(time_of_day = fct_recode(time_of_day, "0" = "AM", "1" = "PM"),
         room_temp = fct_recode(room_temp, "0" = "85", "1" = "90", "2" = "95"),
         group = fct_recode(group, "0" = "Ctrl", "1" = "Trt"),
         loc = fct_recode(loc, "0" = "body_temp_R", "1" = "body_temp_W")) %>%
  mutate(body_temp_diff = body_temp_T - body_temp_P,
         group = fct_relevel(group, c("0", "1")),
         loc = fct_relevel(loc, c("0", "1")),
         num_day = as.numeric(day)
  )

#################################
##### Flank-Rectal Modeling #####
#################################

#model from naively regressing flank vs rectal temperature
simple_fr_mdl <- lm(body_temp_W_adj~body_temp_R_adj, data = chicken_trans)
summary(simple_fr_mdl)

#identifying outliers
count(chicken_trans, body_temp_diff_fac)
#removing outliers
chicken_trans_noout <- filter(chicken_trans, 
                              body_temp_diff_adj <= 0.1, 
                              body_temp_diff_adj >= -0.3)

#random-effects model
rand_fr_mdl <- lme(body_temp_diff_adj~day, random =~1|id
                       , data = chicken_trans_noout)
summary(rand_fr_mdl)
anova(rand_fr_mdl)

##############################
##### Flank-Rectal Plots #####
##############################

#scatterplot for flank and rectal temp
ggplot(chicken_trans, aes(x = body_temp_R_adj, y = body_temp_W_adj)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, col = "red") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  xlab("Rectal Body Temp") + ylab("Flank Body Temp") +
  ggtitle("Flank vs Rectal Body Temp") +
  theme(text=element_text(size=21)) +
  scale_x_continuous(breaks = round(seq(36.5, 42.5, by = 0.5), 1)) +
  theme_bw()

#plot of residuals from simple model
simple_fr_preds <- tibble(fit = predict(simple_fr_mdl),
                          res = residuals(simple_fr_mdl))
ggplot(simple_fr_preds, aes(x = fit, y = res)) +
  geom_point() +
  geom_hline(yintercept = mean(simple_fr_preds$res),
             color = "red", linetype = "dashed") +
  xlab("Fitted Values") + ylab("Residuals") +
  ggtitle("Residuals vs Fitted Values for Flank~Rectal Regression") +
  theme(text=element_text(size=21)) +
  theme_bw()  # Add theme for cleaner look

#histogram of flank-rectal temp
ggplot(chicken_trans, aes(x = body_temp_diff_fac)) +
  geom_bar() +
  xlab("Flank - Rectal Temp") +
  ggtitle("Histogram of Flank-Rectal Temp Difference") +
  theme_bw()

#predicted difference values vs rectal temperatures
chicken_trans_noout$fit <- predict(rand_fr_mdl)
ggplot(chicken_trans_noout, aes(x = body_temp_R_adj, y = fit)) +
  geom_point() +
  geom_hline(yintercept = mean(chicken_trans_noout$fit), color = "red", linetype = "dashed") +
  xlab("Rectal Temperature") + ylab("Fitted Flank-Rectal Difference") +
  ggtitle("Fitted Values for Flank-Rectal Difference from Mixed-Effects Model") +
  theme(text=element_text(size=21)) +
  theme_bw()

######################################
### Cell Body Temperature Modeling ###
######################################

#fixed effects
for(i in 1:2){
  
  loop_loc <- c("R", "F")[i]
  loop_loc_full <- c("Rectal", "Flank")[i]
  
  for(j in 1:4){
  
    #1-way ANOVA
    var1 <- c("group", "time_of_day", "room_temp", "day")[j]
    form1 <- paste0("body_temp_T ~ ", var1)
    one_way_mdl <- lm(as.formula(form1), data = filter(chicken_agg, loc == i-1))
    
    #printing anova table
    cat(paste("\n#############################", 
                "\nANOVA FOR 1-WAY", toupper(var1), 
                toupper(loop_loc_full), "MODEL", 
                "\n#############################\n", sep = " "))
    print(anova(one_way_mdl))   
    
    #assigning 1-way model to variables of the form 
    #VAR1_mdl_LOC for access outside of loop   
    nam <- paste(c(var1, "mdl", loop_loc), collapse = "_")
    assign(nam, one_way_mdl)
    
    for(k in which(1:4 > j)){
      
      #2-way ANOVA
      var2 <- c("group", "time_of_day", "room_temp", "day")[k]
      form2 <- paste0("body_temp_T ~ ", var1, "*", var2)
      two_way_mdl <- lm(as.formula(form2), data = filter(chicken_agg, loc == i-1))     

      #printing anova table
      cat(paste("\n#############################", 
                  "\nANOVA FOR 2-WAY", toupper(var1), "*", 
                  toupper(var2), toupper(loop_loc_full), "MODEL", 
                  "\n#############################\n", sep = " "))
      print(anova(two_way_mdl))   
      
      #assigning 2-way model to variables of the form 
      #VAR1_VAR2_mdl_LOC for access outside of loop
      nam <- paste(c(var1, var2, "mdl", loop_loc), collapse = "_")
      assign(nam, two_way_mdl)
      
    }
  
  } 
}

#fitting random-effect models
chicken_air_mdls <- rep(list(vector("list", 2)), 3)
for(i in 1:3){
  
  loop_temp <- c("85", "90", "95")[i]
  
  for(j in 1:2){
    
    loop_loc <- c("R", "F")[j]
    loop_loc_full <- c("Rectal", "Flank")[j]   
    loop_dat <- filter(chicken_agg, room_temp == i-1, loc == j-1)
    
    chicken_air_mdl <- lme(body_temp_T~rep+day+group+time_of_day,
                         random = ~1|id, 
                         data = loop_dat)   
    
    #printing model summaries
    cat(paste("\n#############################", 
                "\nSUMMARY FOR", loop_temp, "DEGREE", 
                toupper(loop_loc_full), "RANDOM-EFFECTS MODEL", 
                "\n#############################\n", sep = " "))
    print(summary(chicken_air_mdl))
    cat(paste("\n#############################", 
                "\nANOVA FOR", loop_temp, "DEGREE", 
                toupper(loop_loc_full), "RANDOM-EFFECTS MODEL", 
                "\n#############################\n", sep = " "))
    print(anova(chicken_air_mdl))
    
    #assigning model to variable of the form chicken_air_mdl_TEMP_LOC
    #to be accessed outside of loop
    nam <- paste(c("chicken_air_mdl", loop_temp, loop_loc), collapse = "_")
    assign(nam, chicken_air_mdl)
    
    #storing fits for later
    chicken_air_mdls[[i]][[j]] <- chicken_air_mdl
    
  }
  
}

#######################################
##### Cell Body Temperature Plots #####
#######################################

#EDA: plotting body temperature for explanatory variables
#function to get legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#plot to obtain universal legend
legend_plot <- ggplot(chicken_agg, aes(x = group, y = body_temp_T, fill = group)) +
  geom_boxplot(show.legend = T) +
  scale_fill_brewer(type = "qual", direction = -1, 
                    name = "Exp. \nGroup", labels = c("Ctrl", "Trt")) +
  scale_x_discrete(labels = c("0" = "Ctrl", "1" = "Trt")) +
  ggtitle("Experimental Group") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border=element_rect(fill=NA)) +
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank")))

legend <- get_legend(legend_plot)

#cell temp by group
plot1 <- ggplot(chicken_agg, aes(x = group, y = body_temp_T, fill = group)) +
  geom_boxplot(show.legend = F) +
  scale_fill_brewer(type = "qual", direction = -1,
                    name = "Exp. \nGroup", labels = c("Ctrl", "Trt")) +
  scale_x_discrete(labels = c("0" = "Ctrl", "1" = "Trt")) +
  ggtitle("Experimental Group") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border=element_rect(fill=NA)) +
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank")))

#cell temp by time-of-day
plot2 <- ggplot(chicken_agg, aes(x = time_of_day, y = body_temp_T, fill = group)) +
  geom_boxplot(show.legend = F) +
  scale_fill_brewer(type = "qual", direction = -1,
                    labels = c("0" = "Ctrl", "1" = "Trt")) +
  scale_x_discrete(labels = c("0" = "Morning", "1" = "Afternoon")) +
  ggtitle("Time of Day") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border=element_rect(fill=NA)) +
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank")))

#cell temp by room temp
plot3 <- ggplot(chicken_agg, aes(x = room_temp, y = body_temp_T, fill = group)) +
  geom_boxplot(show.legend = F) +
  scale_fill_brewer(type = "qual", direction = -1, 
                    labels = c("0" = "Ctrl", "1" = "Trt")) +
  scale_x_discrete(labels = c("0" = "85", "1" = "90", "2" = "95")) +
  ggtitle("Room Temperature") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border=element_rect(fill=NA)) +
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank")))
  
#cell temp by day
plot4 <- ggplot(chicken_agg, aes(x = day, y = body_temp_T, fill = group)) +
  geom_boxplot(show.legend = F) +
  scale_fill_brewer(type = "qual", direction = -1,
                    labels = c("0" = "Ctrl", "1" = "Trt")) +
  ggtitle("Day") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border=element_rect(fill=NA)) +
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank")))

#cell temp by rep
plot5 <- ggplot(chicken_agg, aes(x = rep, y = body_temp_T, fill = group)) +
  geom_boxplot(show.legend = F) +
  scale_fill_brewer(type = "qual", direction = -1, labels = c("0" = "Ctrl", "1" = "Trt")) +
  ggtitle("Rep") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border=element_rect(fill=NA)) +
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank")))

#arranging plots on common grid
grid.arrange(plot1, plot2, plot3, plot4, plot5, legend,
             top = textGrob("Cell Body Temperature by Factor", gp = gpar(fontsize = 18)),
             left = textGrob("Cell Body Temperature", gp = gpar(fontsize = 18), rot = 90),  
             ncol = 3, nrow = 2)

#analysis: boxplots of fits
#converting fitted values to vector
cell_temp_fits <- chicken_air_mdls %>%
  flatten %>%
  map("fitted") %>%
  map(as_tibble) %>%
  bind_rows %>%
  .$id

#creating dataframe for boxplots
cell_temp_dfs <- lapply(0:2, function(i){
                      lapply(0:1, function(j){
                        return(filter(chicken_agg, room_temp == i, loc == j))
                      })
        }) %>%
  flatten %>%
  bind_rows %>%
  mutate(fit = cell_temp_fits) %>%
  select(loc, room_temp, day, time_of_day, group, fit)

#fitted group effect
ggplot(data = cell_temp_dfs, aes(x = group, y = fit, fill = room_temp)) +
  geom_boxplot() +
  scale_fill_brewer(type = "seq", palette = 8, direction = 1,
                    name = "Room \nTemp.", labels = c("85", "90", "95")) + 
  scale_x_discrete(labels = c("0" = "Ctrl", "1" = "Trt")) + 
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank"))) +
  xlab("Experimental Group") + ylab("Fitted Values") +
  ggtitle("Fitted Group Effect for Rectal and Flank Cell Temperatures") +
  theme_bw()

#fitted day effect
ggplot(data = cell_temp_dfs, aes(x = day, y = fit, fill = room_temp)) +
  geom_boxplot() +
  scale_fill_brewer(type = "seq", palette = 8, direction = 1,
                    name = "Room \nTemp.", labels = c("85", "90", "95")) + 
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank"))) +
  xlab("Day") + ylab("Fitted Values") +
  ggtitle("Fitted Day Effect for Rectal and Flank Cell Temperatures") +
  theme_bw()

#fitted time-of-day effect
ggplot(data = cell_temp_dfs, aes(x = group, y = fit, fill = room_temp)) +
  geom_boxplot() +
  scale_fill_brewer(type = "seq", palette = 8, direction = 1,
                    name = "Room \nTemp.", labels = c("85", "90", "95")) + 
  scale_x_discrete(labels = c("0" = "Morning", "1" = "Afternoon")) + 
  facet_wrap(~loc, labeller = labeller(loc = c(`0` = "Rectal", `1` = "Flank"))) +
  xlab("Time of Day") + ylab("Fitted Values") +
  ggtitle("Fitted Time-of-Day Effect for Rectal and Flank Cell Temperatures") +
  theme_bw()
