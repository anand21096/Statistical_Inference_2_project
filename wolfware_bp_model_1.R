#### Installing required packages if not installed in the system ####
#list.of.packages = c("ggplot2","dpylr","plotly")
#new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

#### 1.1 ####
data_s = read.csv('framingham_data.csv')

#### 1.2 ####
#Getting the classification on number of non-smokers and smokers.
smoker_no = sum(data_s$currentSmoker == 1)
non_smoker_no = sum(data_s$currentSmoker == 0)

#Indexes for non-smoker and smoker. 
smoker_ind = (data_s$currentSmoker == 1)
non_smoker_ind = (data_s$currentSmoker == 0)

########################################################

#Defining dataset for smoker and non smoker.
#Later, we used in mean and variance calculation.
smoker = data_s$sysBP[smoker_ind]
non_smoker = data_s$sysBP[non_smoker_ind]

#Setting alpha value equals to 0.05
alpha = 0.05

#Computing the variance of non-smoker and smoker.
non_smoker_var = var(non_smoker)
smoker_var = var(smoker)

#getting the total degrees of freedom
df = non_smoker_no+smoker_no-2

#Computing the mean of non-smoker and smoker.
non_smoker_mean = mean(non_smoker)
smoker_mean = mean(smoker)


#Getting the pooled variance
p_var = ((smoker_no-1)*smoker_var + (non_smoker_no-1)*non_smoker_var)/(non_smoker_no+smoker_no-2)

#t-stat calculation 
t_stat = (non_smoker_mean - smoker_mean)/(sqrt(p_var)*sqrt(1/smoker_no + 1/non_smoker_no))

#checking whether the value is less than alpha
p_value = 2*pt(t_stat,df,lower.tail = FALSE)
if (p_value > alpha){
  print("Here, we fail to reject the null hypothesis using p-value test while assuming equal variance")
} else{
  print("Here, we reject the null hypothesis using p-value test while assuming equal variance")
}

#We compute standard deviation and t-critical value & use that in confidence interval
p_sd = sqrt(p_var*(1/smoker_no + 1/non_smoker_no))
t_critical = qt(alpha/2,df,lower.tail = FALSE)
upper = (non_smoker_mean - smoker_mean) + t_critical*p_sd
lower = (non_smoker_mean - smoker_mean) - t_critical*p_sd
if (0>lower & 0<upper){
  print("we fail to reject null hypothesis using confidence interval test while assuming equal variance")

} else{
  print("we reject null hypothesis using confidence interval test while assuming equal variance")
}

#Now, we calculate assuming unequal variance (with Satterthwaite approximation for the degrees of freedom)
#Calculating Satterthwaite approximation degree of freedom
df_sw = floor((smoker_var/smoker_no + non_smoker_var/non_smoker_no)**2)/(((smoker_var/smoker_no)**2)/(smoker_no-1) +  ((non_smoker_var/non_smoker_no)**2)/(non_smoker_no-1))

#We first calculate the t value and use it in p-value test
t_stat = (non_smoker_mean - smoker_mean)/(sqrt(smoker_var/smoker_no + non_smoker_var/non_smoker_no))
p_value = 2*pt(t_stat,df_sw,lower.tail = FALSE)
if (p_value < alpha){
  print("we reject the null hypothesis using p-value test while assuming unequal variance")
} else{
  print("we fail to reject the null hypothesis using p-value test while assuming unequal variance")
}

#First, we calculate t-critical and standard deviation.
#Later, We use it in upper bound and lower bound calculation
t_critical = qt(alpha/2,df_sw,lower.tail = FALSE)
p_sd = sqrt(smoker_var/smoker_no + non_smoker_var/non_smoker_no)
upper = (non_smoker_mean - smoker_mean) + t_critical*p_sd
lower = (non_smoker_mean - smoker_mean) - t_critical*p_sd
if (0>lower & 0<upper){
  print("we fail to reject null hypothesis using confidence interval test while assuming unequal variance")
} else{
  print("we reject null hypothesis using confidence interval test while assuming unequal variance")
}


#### 1.3 ####
#QQ-plot for non-smoker
qqnorm(data_s$sysBP[non_smoker_ind], pch = 10, col = "blue", frame = FALSE, main="QQ plot for systolic BP")
qqline(data_s$sysBP[non_smoker_ind],lwd = 4, col = "maroon")

#QQ-plot for smoker
qqnorm(data_s$sysBP[smoker_ind], pch = 10, col = "blue", frame = FALSE, main="QQ plot for systolic BP")
qqline(data_s$sysBP[smoker_ind],lwd = 4, col = "maroon")

#### 1.4 ####
#From QQ-plot, we conclude that the data is not normal. 
#Hence, we use Mann-Whitney test. 
#First, we pool the data in the order of smoker and non smoker. 
#As smoker data is minority, we rank the smoker data and used it in wilcox test
df = c(smoker, non_smoker)
rank_smoker = rank(df)[1:smoker_no]
pwilcox(sum(rank_smoker), smoker_no, non_smoker_no)
#Since this value is above 0.05, we do not reject the null hypothesis. This is different from 
# the previous tests assuming normality.

