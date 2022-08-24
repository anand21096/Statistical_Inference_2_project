#### Installing required packages if not installed in the system ####
#list.of.packages = c("ggplot2","dpylr","plotly")
#new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

#Importing required libraries
library(ggplot2)
library(dplyr)
library(plotly)


#defining equal variance test function
Variance_eq_Testing <- function(a,b){
  #Threshold level
  threshold_alpha = 0.05
  a_len = length(a) 
  b_len = length(b)
  a_variance = var(a)
  b_variance = var(b)
  mean_a = mean(a)
  mean_b = mean(b)
  
  #getting degrees of freedom 
  degrees_of_freedom = a_len+b_len-2
  #pooled variance
  pool_var = ((a_len-1)*a_variance + (b_len-1)*b_variance)/(degrees_of_freedom)
  #defining t-statistic
  t_stat = (mean_a-mean_b)/(sqrt(pool_var)*sqrt(1/a_len+1/b_len))
  #function returns TRUE if null hypothesis is rejected & FALSE otherwise
  return(abs(t_stat)>qt(threshold_alpha/2, degrees_of_freedom,  lower.tail = FALSE))
}

#defining unequal variance test function
Variance_uneq_Testing <- function(p,q){
  threshold_alpha = 0.05
  mean_p = mean(p)
  mean_q = mean(q)
  p_len = length(p) 
  q_len = length(q)
  p_variance = var(p)
  q_variance = var(q)
  
  #calculating degrees of freedom
  degrees_of_freedom = ((((p_variance/p_len) + (q_variance/q_len))^2) / (((p_variance/p_len)^2/(p_len-1)) + (q_variance/q_len)^2/(q_len-1)))
  
  #getting t-stat value
  t_stat = (mean_p - mean_q) / sqrt(p_variance/p_len + q_variance/q_len)
  
  #this function returns TRUE if null hypothesis is rejected. If not, it returns FALSE
  return(abs(t_stat) > qt(threshold_alpha/2, degrees_of_freedom, lower.tail = FALSE))
}


#simulate 10000 instances and perform t-test after defining simulation function.
simulation <- function(sample_size1, sample_size2, mean, std_dev, t_test_value){
  #taking 10000 simulations
  iterations = 10000
  #list for storing status of each simulation(accept or reject)
  list_res = c()
  #iterating over the 1000 instances
  for (i in 1:iterations){
    #simulating a sample from normal distribution with mean difference "mean" and standard deviation "std_dev"
    sample1 = rnorm(sample_size1, mean, std_dev)
    #simulating a sample from normal distribution with mean 0 and standard deviation 1
    sample2 = rnorm(sample_size2, 0, 1)
    #selecting test according to input by user (equal variance vs unequal variance)
    if (t_test_value == "equal"){
      #returning list holding TRUE or FALSE values
      list_res[i] = Variance_eq_Testing(sample1, sample2)
    } else {
      #returning list holding TRUE or FALSE values 
      list_res[i] = Variance_uneq_Testing(sample1, sample2)
    }
  }
  #rejection of null hypothesis is indicated by true and FALSE indicates otherwise
  result = length(list_res[list_res == TRUE])/iterations
  #returns the fraction of TRUE values
  return(result)
}

#initializing dataframe to hold the return of simulation function for each combination
output <- as.data.frame(matrix(0,1,6))

difference_in_mean <- c(0,-5,5,-1,1)
std_dev <- c(1,2,3)
n <- c(10,30,70)

#naming the columns of the output dataframe accordingly
names(output) <- c("type_t_test", "sample_size", "std_dev", "Mean_Diff","Threshold_type","threshold_value")

#iterating over all possible permutations and combinations
for (p in std_dev){
  for (q in n){
    for (r in n){
      for (s in difference_in_mean){
        std_dev_val <- paste("[",p,",",1,"]")
        size <- paste("[",q,",",r,"]")
        #performing equal variance test & storing data in dataframe
        output[nrow(output)+1,] <- c("Equal Variance",size, std_dev_val, s, ifelse(s == 0, "alpha", "power"), simulation(q,r,s,p,"equal"))
        #performing unequal variance test & storing data in dataframe
        output[nrow(output)+1,] <- c("Unequal Variance", size, std_dev_val, s, ifelse(s == 0, "alpha", "power"), simulation(q,r,s,p,"unequal"))
      }
    }
  }
}
#displaying the final output
##removing the initial dummy row
output <- output[-c(1),]
head(output, n = 25)

#dividing data based on standard deviation of samples
data_one <- output[output$std_dev == "[ 1 , 1 ]",]
data_two <- output[output$std_dev == "[ 2 , 1 ]",]
data_three <- output[output$std_dev == "[ 3 , 1 ]",]

#plotting graph with the data_one, data_two & data_three
graph_1 <- ggplot(data.frame(data_one), aes(x=sample_size,y=threshold_value))+geom_point(aes(shape = Threshold_type, color = Mean_Diff)) + xlab("(n1, n2)") + ylab("Type I or Type II Error") + ggtitle("std_dev_1 = 1")
graph_1 + facet_wrap(~type_t_test)
graph_2 <- ggplot(data.frame(data_two), aes(x=sample_size,y=threshold_value))+geom_point(aes(shape = Threshold_type, color = Mean_Diff)) + xlab("(n1, n2)") + ylab("Type I or Type II Error") + ggtitle("std_dev_1 = 2")
graph_2 + facet_wrap(~type_t_test)
graph_3 <- ggplot(data.frame(data_three), aes(x=sample_size,y=threshold_value))+geom_point(aes(shape = Threshold_type, color = Mean_Diff)) + xlab("(n1, n2)") + ylab("Type I or Type II Error") + ggtitle("std_dev_1 = 3")
graph_3 + facet_wrap(~type_t_test)

#The End