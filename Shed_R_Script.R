# setwd('~/')

# install.packages("pracma")
library(tidyverse)
library(pracma)

### Shedding data
shed_data <- read_csv('FS19_ForAUCAnalysis.csv')


sum_shed <- shed_data %>%
  group_by(Animal) %>%
  summarise(AUC=trapz(Day, Shed),
            AULC=trapz(Day, log(Shed+1)),
            sum=sum(Shed), 
            maxshed=Shed[which.max(Shed)], 
            day_max=Day[which.max(Shed)], 
            pos_samples=sum(Shed > 0), 
            StudyVar=unique(StudyVar), 
            Strain=unique(Strain))



sum_shed %>%
  ggplot(aes(x=Strain, y=AULC), fill=StudyVar) +
  geom_boxplot()+
  geom_jitter(width= 0.22, size=3, aes(color=StudyVar, fill=StudyVar))+
  ggtitle("Fecal Shedding AULC by Strain") + xlab("Strain") + labs(color="Trial", fill="Trial")




summary(aov(data = sum_shed, formula = AULC~Strain))
aovDat <- aov(data = sum_shed, formula = AULC~Strain)
TukeyHSD(aovDat)



# jules add
shed_time_sum <- 
  shed_data %>% 
  group_by(Day, Strain) %>% 
  summarise(avg_shed=mean(Shed), 
            avg_log_shed=mean(log(Shed + 1)))


shed_time_sum %>% 
  ggplot(aes(x=Day, y=avg_log_shed, group=Strain, color=Strain)) + 
  geom_point() +
  geom_line()



#### SWABS
swab_data <- read_csv('RamsSwabs.csv')


sum_swab <- swab_data %>% group_by(Animal) %>%
  summarise(AUC=trapz(Day, Swab),
            AULC=trapz(Day, log(Swab+1)),
            sum=sum(Swab), 
            maxswab=Swab[which.max(Swab)], 
            day_max=Day[which.max(Swab)], 
            pos_samples=sum(Swab > 0), 
            StudyVar=unique(StudyVar), 
            Strain=unique(Strain))


sum_swab %>%
  ggplot(aes(x=Strain, y=AULC), fill=StudyVar) +
  geom_boxplot()+
  geom_jitter(width= 0.22, size=3, aes(color=StudyVar, fill=StudyVar))+
  ggtitle("Fecal Swab AULC by Strain") +
  xlab("Strain") + 
  labs(color="Trial", fill="Trial")




summary(aov(data = sum_swab, formula = AULC~Strain))
aovDat <- aov(data = sum_swab, formula = AULC~Strain)
TukeyHSD(aovDat)



# Jules Add

swab_data %>%
  ggplot(aes(x=Day, y=log(Swab + 1), group=Animal, color=Strain)) + geom_line()

swab_time_sum <- 
  swab_data %>%
  group_by(Day, Strain) %>% 
  summarise(avg_swab=mean(Swab), 
            avg_log_swab=mean(log(Swab + 1))) 


swab_time_sum %>% 
  ggplot(aes(x=Day, y=avg_log_swab, group=Strain, color=Strain)) + 
  geom_point() +
  geom_line()

shed_time_sum



swab_data <- swab_data %>%
  mutate(type='swab',
         CFUs=Swab, 
         Strain=sub('(.*)W','\\1',Strain), 
         Strain=sub('T14588', 'TW14588', Strain)) %>% 
  select(-Swab)

shed_data <- shed_data %>%
  mutate(type='shed',
         CFUs=Shed) %>% 
  select(-Shed)

all_dat <- bind_rows(shed_data, swab_data)

test <- all_dat %>% 
  pivot_wider(names_from = type, values_from=CFUs) %>% 
  mutate(site_dif=shed-swab, 
         log_site_diff=log(shed + 1) - log(swab + 1))

test %>%
  ggplot(aes(x=Day, y=log_site_diff, group=Animal)) +
  geom_point() + geom_line() + facet_wrap(~Strain)

test_sum <- test %>% group_by(Day, Strain) %>% 
  summarise(avg_log_site_diff = mean(log_site_diff), 
            stderr=std_err(log_site_diff), 
            .groups='drop')

test_sum %>%
  ggplot(aes(x=Day, y=avg_log_site_diff, group=Strain, color=Strain))+ 
  geom_point()+
  geom_line() + 
  geom_errorbar(aes(ymin=avg_log_site_diff - stderr, 
                    ymax=avg_log_site_diff + stderr), width=.2) + 
  facet_wrap(~Strain)



min(test$site_dif)


all_sum <- all_dat %>% group_by(Day, Strain, type) %>% 
  summarise(avg_log_CFU=mean(log(CFUs + 1)), 
            stderr=std_err(log(CFUs + 1)), 
            .groups='drop')


all_sum %>% ggplot(aes(x=Day,
                       y=avg_log_CFU,
                       group=interaction(Strain, type), 
                       color=type)) + 
  geom_errorbar(aes(ymin=avg_log_CFU - stderr, 
                    ymax=avg_log_CFU + stderr),
                width=.2,
                alpha=.7)+
  geom_point() + geom_line() + facet_wrap(~Strain, nrow = 1) + 
  ylim(0,12)

