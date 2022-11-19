
#Code for submission

#Install libraries
library(tidyverse)
library(deSolve)
library(abind)
library(pracma)
library(ggpubr)
library(plotly)
library(ggsci)

set.seed(1004)

#set fixed parameter values
max_dur_years <- (0.25 + hill_down_function(10^2) + 0.75) #max duration in yrs 
max_dur <- max_dur_years*365
num_strains_A=150
s_A=10^-4
num_strains_B=35
s_B=4*10^-3
num_strains_C=10
s_C=10^-1

s =  logseq(10^-1,10^-5, n=50)
x = logseq(0.5,0.012, n=50)

vl_A = 7-(0:(num_strains_A))*x[(which.min(abs(s-s_A)))] #spvl approximatley goes between 
vl_B = 7-(0:(num_strains_B))*x[(which.min(abs(s-s_B)))]
vl_C = 7-(0:(num_strains_C))*x[(which.min(abs(s-s_C)))]



fullmodelA <- quasispecies_model_full(m=num_strains_A,M=num_strains_A+1,s=s_A, max_dur=20000, dt=1, mu=mu)
fullmodelB <- quasispecies_model_full(m=num_strains_B,M=num_strains_B+1,s=s_B, max_dur=max_dur, dt=1, mu=mu)
fullmodelC <- quasispecies_model_full(m=num_strains_C,M=num_strains_C+1,s=s_C, max_dur=max_dur, dt=1, mu=mu)

equilA <- cbind(quasispecies_model_equilibrium_full(m=num_strains_A,M=num_strains_A+1,s=s_A, mu=mu), spvl=vl_A, mod=str(num_strains_A-1))
equilB <- cbind(quasispecies_model_equilibrium_full(m=num_strains_B,M=num_strains_B+1,s=s_B, mu=mu), spvl=vl_B, mod=str(num_strains_B-1))
equilC <- cbind(quasispecies_model_equilibrium_full(m=num_strains_C,M=num_strains_C+1,s=s_C, mu=mu), spvl=vl_C, mod=str(num_strains_C-1))
equil <- rbind(equilA, equilB, equilC) #dist at equilibrium

vl_time_A <- viralload_over_time(num_strains_A, fullmodelA, vl_A)
vl_time_B <- viralload_over_time(num_strains_B, fullmodelB, vl_B)
vl_time_C <- viralload_over_time(num_strains_C, fullmodelC, vl_C)

#extended version of the model, 30,000 days: 



FOI_A = force_of_infection_calculation(time=max_dur, spvl=vl_A, muts_freq=fullmodelA, E=0, num_strains_A+1 )
FOI_B = force_of_infection_calculation(time=max_dur, spvl=vl_B, muts_freq=fullmodelB, E=0, num_strains_B+1)
FOI_C = force_of_infection_calculation(time=max_dur, spvl=vl_C, muts_freq=fullmodelC, E=0, num_strains_C+1)

durations_A <- FOI_A[[2]]
durations_B <- FOI_B[[2]]
durations_C <- FOI_C[[2]]

spvl_overTimeA = FOI_A[[3]]
spvl_overTimeB = FOI_B[[3]]
spvl_overTimeC = FOI_C[[3]]

spvl_A = c()
spvl_B = c()
spvl_C = c()

for (i in 1:(num_strains_A+1)){
  d = durations_A[i]*365
  spvl_A[i] = mean(spvl_overTimeA[1:d,i])
}

for (i in 1:(num_strains_B+1)){
  d = durations_B[i]*365
  spvl_B[i] = mean(spvl_overTimeB[1:d,i])
}

for (i in 1:(num_strains_C+1)){
  d = durations_C[i]*365
  spvl_C[i] = mean(spvl_overTimeC[1:d,i])
}


#between host model at equilibria: 
equil_prev_A <- bh_equilbria(FOI_A, durations_A, dt=1, spvl_A,num_strains_A+1, "A")
equil_prev_B <- bh_equilbria(FOI_B, durations_B, dt=1, spvl_B,num_strains_B+1, "B")
equil_prev_C <- bh_equilbria(FOI_C, durations_C, dt=1, spvl_C,num_strains_C+1, "C") #R0 too small
equil_prev <- rbind(equil_prev_A[[1]], equil_prev_B[[1]], equil_prev_C[[1]])


#between host model 
time_taken = 50 #time step
dt=time_taken/365
AICint = max_dur_years  #Duration of initial conditions: the max duration of the infection
iT0 = ceil(AICint/dt) #Time in days
num_host_types = 1
inds <- c(TRUE,rep(FALSE, time_taken-1))
fA = FOI_A[[1]][,,inds]
fB = FOI_B[[1]][,,inds]
fC = FOI_C[[1]][,,inds]
modA = list()
modB = list()
modC = list()

# #start with high viral load/average viral load/low viral load
initialVL_A = data.frame(sample(1:50, 10), initialVL_A = sample(50:90, 10), initialVL_A = sample(110:150, 10))
initialVL_B= data.frame(sample(1:10, 10), initialVL_A = sample(11:23, 10), initialVL_A = sample(26:36, 10))
initialVL_C = data.frame(sample(1:3, 3), initialVL_A = sample(4:7, 3), initialVL_A = sample(8:10, 3))

for (i in 1:3){
 modA[[i]] = between_host_model(f=fA,
                                 tmax=150,
                                 tstep=time_taken,
                                 max_dur=max_dur,
                                 initCons=initialVL_A[,i],
                                 num_strains=num_strains_A+1,
                                 num_hosts=1,
                                 durations_i=durations_A,
                                 host_dist=1,
                                 transmissibility=1,
                                 start_size=1)
}

for (i in 1:3){
  modB[[i]] = between_host_model(f=fB,
                                 tmax=100,
                                 tstep=time_taken,
                                 max_dur=max_dur,
                                 initCons=initialVL_B[,i],
                                 num_strains=num_strains_B+1,
                                 num_hosts=1,
                                 durations_i=durations_B,
                                 host_dist=1,
                                 transmissibility=1,
                                 start_size=1)
}

for (i in 1:3){
  modC[[i]] = between_host_model(f=fC,
                                 tmax=100,
                                 tstep=time_taken,
                                 max_dur=max_dur,
                                 initCons=initialVL_C[,i],
                                 num_strains=num_strains_C+1,
                                 num_hosts=1,
                                 durations_i=durations_C,
                                 host_dist=1,
                                 transmissibility=1,
                                 start_size=3)
}


#Alternatively load in the data: 
load("~/Documents/DPhilWork/Data/MSB work/Results_MSBmodel/BH_numerical_results_oct/model_bh_A.RData")
load("~/Documents/DPhilWork/Data/MSB work/Results_MSBmodel/BH_numerical_results_oct/model_bh_B.RData")
load("~/Documents/DPhilWork/Data/MSB work/Results_MSBmodel//BH_numerical_results_oct/model_bh_C.RData")

modA1 = modA[[1]][[3]]
modA2 = modA[[2]][[3]]
modA3 = modA[[3]][[3]]
modB1 = modB[[1]][[3]]
modB2 = modB[[2]][[3]]
modB3 = modB[[3]][[3]]
modC1 = modC[[1]][[3]]
modC2 = modC[[2]][[3]]
modC3 = modC[[3]][[3]]

host_effect <- analytical_sol_heteropop(50, num_strains_A, scenarioExtra_fullmodelA,  s_x=s_A, spvl_A,E, "1log", dist)
host_effect_df <- host_effect[[4]]

#Heritability 
herit_A_m1 <- heritability_m1(500,abs(equil_prev_A[[1]]$prevalence), spvl_A, num_strains_A+1, FOI_A[[4]])
herit_A_m2 <- heritability_m2(500,abs(equil_prev_A[[1]]$prevalence), spvl_A, num_strains_A+1, FOI_A[[4]])


herit_50_A_halflog_m1 <- heritability_with_host_effect_m1(P=500, num_strains=num_strains_A, 50, "halflog", "A")
herit_50_A_halflog_m2 <- heritability_with_host_effect_m2(P=500, num_strains=num_strains_A, 50, "halflog", "A")
herit_noeffect_A_m1 <- data.frame(mean=herit_A_m1[[1]], std=herit_A_m1[[2]], num_strains=num_strains_A,
                                  host_types=0, host_size="Anone", model="A")
herit_noeffect_A_m2 <- data.frame(mean=herit_A_m2[[1]], std=herit_A_m2[[2]], num_strains=num_strains_A,
                                  host_types=0, host_size="Anone", model="A")

#Combined dataframe. 
heritability_df <- rbind(herit_noeffect_A_m1 , herit_noeffect_A_m2, 
                         herit_50_A_halflog_m1[[1]], herit_50_A_halflog_m2[[1]])
heritability_df$methods <- c("Method 1", "Method 2", "Method 1", "Method 2")
heritability_df$host_size <- factor(c("Homogenous", "Homogenous", "Heterogeneous", "Heterogeneous"),
                                    levels = c("Homogenous","Heterogeneous"))


######### Main paper figures #########


#figure 1: 
fig1a <- data.frame(fitness_cost=s, viral_cost=x) %>% ggplot() + aes(x=fitness_cost, y=viral_cost) + 
  theme_classic() + geom_line(size=1.2) + scale_y_log10(limit=c(0.01,0.5), 
                                                        breaks=c(0.01, 0.1, 0.5)) +
  xlab("Fitness cost") + ylab("Viral load cost") + theme(axis.text = element_text(size=15),
                                                         axis.title = element_text(size=15))+
  scale_x_log10(breaks=c(10^-1, 10^-2, 10^-3, 10^-4, 10^-5) ,labels=c(expression(10^{-1}),expression(10^{-2}),
                                                                      expression(10^{-3}), expression(10^{-4}),
                                                                      expression(10^{-5}))) + 
  geom_point(aes(y=x[(which.min(abs(s-s_A)))],x=s_A), size=10, colour="#BC3C29FF") +
  geom_point(aes(y=x[(which.min(abs(s-s_B)))],x=s_B), size=10, colour="#0072B5FF") +
  geom_point(aes(y=x[(which.min(abs(s-s_C)))],x=s_C), size=10, colour="#3CB371") +
  annotate("text",label="150 sites", x=5*10^-4, y=x[(which.min(abs(s-s_A)))], size=6, colour="#BC3C29FF") +
  annotate("text",label="35 sites", x=20*10^-3, y=x[(which.min(abs(s-s_B)))], size=6, colour="#0072B5FF") +  
  annotate("text",label="10 sites", x=2*10^-2, y=x[(which.min(abs(s-s_C)))], size=6, colour="#3CB371")   

purple_colours <- c("#DDA0DD", "#BA55D3","#8A2BE2", "#4B0082")
fig1b <- ggplot(df, aes(x=m, group=s, colour=factor(s), y=freq)) + geom_line(size=2) + 
  theme_classic() + xlim(c(0,75)) +   
  scale_color_manual(values=purple_colours,labels=c(expression(10^{-5}),expression(10^{-4}),
                                                    expression(10^{-3}), expression(10^{-2})
  )) + labs(colour="Fitness cost") + 
  xlab("Number of deleterious mutations") + ylab("Density") + theme(legend.position = c(0.4,0.8),
                                                                    axis.text=element_text(size=15),
                                                                    axis.title=element_text(size=15),
                                                                    legend.text=element_text(size=15),
                                                                    legend.title = element_text(size=15)) 


#figure 2:
colour_palette <- c("#BC3C29FF","#0072B5FF","#3CB371")
equil$mod <- factor(equil$mod,levels = c("150","35","10"))

fig2a <- ggplot(equil, aes(x=mutations, y=equilibrium, group=mod))  + xlim(c(0,60))+
  ylab("Density") + theme_classic() + xlab("Virus type") + geom_ribbon(aes(fill=mod, ymin=0, ymax=equilibrium), alpha=1) + geom_line() +
  guides(fill=guide_legend(title="Maximum number of mutations"), colour="none") + theme(axis.text = element_text(size = 15),
                                                                                        axis.title = element_text(size=15),
                                                                                        legend.text = element_text(size=20),
                                                                                        legend.title = element_text(size=20),
                                                                                        legend.position=c(0.58,0.8)) +
  scale_fill_manual(values=colour_palette, breaks=c("10","35","150")) +  scale_colour_manual(values=colour_palette, breaks=c("10","35","150"))

equil_parameterA = quasispecies_model_equilibrium_full(m=num_strains_A,M=num_strains_A+1,s=s_A, mu=mu)
equil_spvl_A = sum(equil_parameterA$equilibrium * spvl_A)   
equil_parameterB = quasispecies_model_equilibrium_full(m=num_strains_B,M=num_strains_B+1,s=5*10^-3, mu=mu)
equil_spvl_B = sum(equil_parameterB$equilibrium * spvl_B)
equil_parameterC = quasispecies_model_equilibrium_full(m=num_strains_C,M=num_strains_C+1,s=10^-1, mu=mu)
equil_spvl_C = sum(equil_parameterC$equilibrium * spvl_C)

fig2c  <- ggplot(spvl_time_A %>% drop_na(), aes(x=time, y=spvl, group=start)) + geom_line(colour="#BC3C29FF") +
  guides(colour="none") + theme_minimal() + geom_hline(yintercept=equil_spvl_A) +
  ylab( expression(paste("Viral load (",log[10],")"))) + ylim(c(2,7)) +
  xlab("Days into  chronic infection") + theme(axis.text=element_text(size=15), axis.title=element_text(size=15))

fig2c  <- ggplot(spvl_time_B, aes(x=time, y=spvl, group=start)) + geom_line(colour="#0072B5FF") +
  guides(colour="none") + theme_minimal() + geom_hline(yintercept=equil_spvl_B) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) + xlim(c(0,3000)) +
  xlab("Days into chronic infection") + ylab( expression(paste("Viral load (",log[10],")")))

fig2b  <- ggplot(spvl_time_C, aes(x=time, y=spvl, group=start)) + geom_line(colour="#3CB371") +
  guides(colour="none") + theme_minimal() + xlim(c(0,300)) + geom_hline(yintercept=equil_spvl_C) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15))+ xlab("Days into chronic infection") +
  ylab( expression(paste("Viral load (",log[10],")")))

ggarrange(fig2a,NULL, fig2b, fig2c, NULL, fig2d,
          ncol=3, nrow=2, widths=c(1,0.1,1,1,0.1,1),
          labels=c("A","","B","C","","D"), font.label = list(size=15))


#figure 3: 
#heat map for NGM values. 
KM_A = FOI_A[[6]]
plot3a = data.frame(tp=colSums(KM_A), starting_mutations = 0:num_strains_A) %>% 
  ggplot() + aes(x=starting_mutations, y=tp) + geom_line(colour=colour_palette[1], size=2) + ylim(c(0,2)) +
  theme_classic() + ylab("Transmission\npotential") + xlab('Initial number of mutations') + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15)) + geom_hline(yintercept = 1, alpha=0.8, linetype=2, colour=colour_palette[1])

KM_B = FOI_B[[6]]
plot3b = data.frame(tp=colSums(KM_B), starting_mutations = 0:num_strains_B) %>% 
  ggplot() + aes(x=starting_mutations, y=tp) +  geom_line(colour=colour_palette[2], size=2) + ylim(c(0,2)) +
  theme_classic() + ylab("Transmission\npotential") + xlab('Initial number of mutations') + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15))+ geom_hline(yintercept = 1, alpha=0.8, linetype=2, colour=colour_palette[2])


KM_C = FOI_C[[6]]
plot3c = data.frame(tp=colSums(KM_C), starting_mutations = 0:num_strains_C) %>% 
  ggplot() + aes(x=starting_mutations, y=tp) + geom_line(colour=colour_palette[3], size=2) + ylim(c(0,2)) +
  theme_classic() + ylab("Transmission\npotential") + xlab('Initial number of mutations') + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15)) + scale_x_continuous(breaks=c(0,5,10))+ geom_hline(yintercept = 1, alpha=0.8, linetype=2, colour=colour_palette[3])


plot3abc <- ggarrange(plot3a, plot3b, plot3c, ncol=1, labels=c("A","B","C"), font.label =list(size=20))



r0_df <- data.frame(r0=c(as.numeric(equil_prev_A[[4]]), 
                         as.numeric(equil_prev_B[[4]]),
                         as.numeric(equil_prev_C[[4]])),
                    mod=c("150 segregating sites", "35 segregating sites", "10 segregating sites"))

colour_palette <- c("#3CB371","#0072B5FF","#BC3C29FF")
r0_df$mod <- factor(r0_df$mod,levels = c("10 segregating sites", "35 segregating sites", "150 segregating sites"))
fig3d <- r0_df %>% ggplot() + aes(x=mod, y=r0, fill=mod) + geom_bar(stat="identity") +
  scale_fill_manual(values=colour_palette) + ylab("Basic reproduction number")+ xlab("") +
  theme_classic() + theme(axis.text.x = element_blank(), legend.title=element_blank(),
                          axis.text.y = element_text(size=15), axis.title.y=element_text(size=15), 
                          legend.position = c(0.3,0.9), 
                          legend.text = element_text(size=15)) 

ggarrange(fig3abs, NULL, figd, ncol=3, widths=c(1,0.2, 1),labels=c("","","D"), font.label=list(size=20))

#figure 4:
bh_colour_palette= c("#1B9E77", "#E7298A","#E6AB02")
fig4d <- equil_prev_A[[1]] %>% ggplot() + aes(x=spvl, y=prevalence) + geom_histogram(stat="identity", fill="#BC3C29FF") + 
  theme_classic() +   ylab("Infection prevalence")+ xlim(c(3,6)) +
  xlab(expression(paste("spVL (",log[10],")"))) +   theme(axis.text = element_text(size = 15),axis.title = element_text(size=15),
                                                          legend.text = element_text(size=15),
                                                          legend.title = element_text(size=15))

modA1 = modA[[1]][[3]]
modA2 = modA[[2]][[3]]
modA3 = modA[[3]][[3]]

dist_high <- data.frame(p=modA[[1]][[3]][,494]*modA[[1]][[1]][494,2], y=spvl_A)
dist_int <- data.frame(p=modA[[2]][[3]][,227]*modA[[2]][[1]][227,2], y=spvl_A)
dist_low <- data.frame(p=modA[[3]][[3]][,291]*modA[[3]][[1]][291,2], y=spvl_A)

fig4a <- ggplot(dist_high, aes(x=y, y=p)) + geom_histogram(stat="identity", fill=bh_colour_palette[1])+
  theme_classic() + xlim(c(3,6)) + xlab(expression(paste("spVL (",log[10],")"))) + ylab("Infection\nprevalence") +
  ylim(c(0,15)) + theme(axis.text = element_text(size = 15),axis.title = element_text(size=15),
                        legend.text = element_text(size=15),
                        legend.title = element_text(size=15)) + ggtitle("45 years into epidemic")

fig4b <- ggplot(dist_int, aes(x=y, y=p)) + geom_histogram(stat="identity", fill=bh_colour_palette[2])+
  theme_classic() + xlim(c(3,6)) + xlab(expression(paste("spVL (",log[10],")"))) + ylab("Infection\nprevalence")+
  ylim(c(0,15))+ theme(axis.text = element_text(size = 15),axis.title = element_text(size=15),
                       legend.text = element_text(size=15),
                       legend.title = element_text(size=15)) + ggtitle("10 years into epidemic")

fig4c <- ggplot(dist_low, aes(x=y, y=p)) + geom_histogram(stat="identity", fill=bh_colour_palette[3])+
  theme_classic() + xlim(c(3,6)) + xlab(expression(paste("spVL (",log[10],")"))) + ylab("Infection\nprevalence")+
  ylim(c(0,15))+ theme(axis.text = element_text(size = 15),axis.title = element_text(size=15),
                       legend.text = element_text(size=15),
                       legend.title = element_text(size=15)) + ggtitle("20 years into epidemic") 

#empty legend plot: 
legend=data.frame(InitialVL=c("High", "Intermediate","Low"), x=1, y=1)
p3 <-ggplot(legend, aes(x=x, y=y, colour=InitialVL)) + geom_point() + 
  scale_color_manual(values=bh_colour_palette, name="Initial VL") +
  lims(x = c(0,0), y = c(0,0)) +
  theme_void()+ 
  theme(legend.position = "top",
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =  15),
        legend.title = element_text(size = 15, face = "bold"))

fig4abs <- ggarrange(p3,fig4a, fig4b, fig4c, ncol = 1, labels=c("","A", "B", "C"),
                      heights=c(0.15,1,1,1))


fig4e =  %>% filter(spvl>3, spvl<6) %>%
  mutate(group=cut(spvl,150),
         vl = recode(group,
                     "(3,3.02]"="3",
                     "(4,4.02]"="4",
                     "(4.98,5]"="5",
                     "(5.98,6]" = "6")) %>%
  group_by(vl) %>% summarise(freq=sum(prevalence)) %>% ggplot(aes(x=vl, y=freq)) +
  geom_bar(stat="identity", fill="#BC3C29FF") + theme_classic() +
  ylab("Infection prevalence")+ xlab("") +
  scale_x_discrete(breaks=c("3","4","5",
                            "6"))  + theme(axis.text = element_text(size = 15),axis.title = element_text(size=15),
                                           legend.text = element_text(size=15),
                                           legend.title = element_text(size=15)) + 
  xlab(expression(paste("spVL (",log[10],")")))

fig4de <- ggarrange(fig4d,fig4e,common.legend=T, labels=c("C","D"), nrow=2)          

ggarrange(fig4abc, fig4de, ncol=2)


#figure 5:



#figure 6:
fig6 <- ggplot(heritability_df, aes(group=methods, x=host_size, y=mean, colour=methods))+
  geom_point(position=position_dodge(width=0.2), size=4) + theme_classic()+
  ylim(c(0,1)) + ylab("Heritability estimate") + xlab("Host population type")+
  labs(colour="Source sampling method") +  geom_errorbar(aes(ymax=mean+std, ymin=mean-std),
                                                         position=position_dodge(width=0.2), width=0.075) +
  theme(axis.text = element_text(size = 20),axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20), legend.position=c(0.3,0.3)) 

#supp figure: 








