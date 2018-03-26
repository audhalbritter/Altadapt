make_data <- function(np, d, ind, sd_ind, sd_pop, sd_hab, sd_dem, SA, constant)
{
  
  #PARAMETER DEFINITIONS:
  
  #np is the number of populations tested
  #d indicates design, from d = 2 (twice as many allopatric as sympatric transplant), to np - 1 (full factorial design)
  #ind is the number of individuals tested by transplant
  
  #sd_ind is the individual - level standard deviation of error
  #sd_pop is the population - level standard deviation of error
  #sd_hab is the standard deviation of habitat effects
  #sd_dem is the standard deviation of deme quality effects
  
  #SA is the sympatric vs. allopatric contrast (i.e., local adaptation)
  #constant is the baseline fitness
  
  #habitat and deme quality effects are drawn in a normal distribution with sd as defined above:
  habeffects=rnorm(np,m = 0,sd = sd_hab)
  demeffects=rnorm(np,m = 0,sd = sd_dem)
  
  #now generate the transplant matrix data:
  fitness = c(); habs = c(); dems = c(); symp = c(); 
  
  for(i in 1:np){ #loop on deme of destination
    for(j in 1:np){ #loop on deme of origin
      
      #draw the population-level error and add the "SA" contrast to sympatric transplants:
      if(i == j){error_pop = rnorm(1,m = SA,sd = sd_pop)}
      if(i != j){error_pop = rnorm(1,m = 0,sd = sd_pop)}
      
      for(k in 1:ind){
        #draw error for this individual:
        error_ind = rnorm(1,m = 0,sd = sd_ind);
        
        #Generate fitness data for sympatric transplant:
        if(i == j){
          fitness = append(fitness,constant + habeffects[i] + demeffects[j] + error_pop + error_ind)
          symp = append(symp,1)
          habs=append(habs,i)
          dems=append(dems,j)					   
        }
        
        #Generate fitness data for allopatric transplant, depending on the experimental design d
        if((i-j >= 1 - np & i-j <= d - np) || (i-j >= 1 & i-j <= d)){
          fitness=append(fitness,constant + habeffects[i] + demeffects[j] + error_pop + error_ind)
          symp = append(symp,0)
          habs=append(habs,i)
          dems=append(dems,j)	
        }
      }  
    }
  }
  #return the transplant experiment data
  return(data.frame(habs, dems, symp, fitness))
}






#we use the same function but we assume 1 individual is sampled and sd_ind is 0 (equivalent to assuming one measure per population)
np = 10	#np is the number of populations tested
d = 2	#d indicates design, from d = 2 (twice as many allopatric as sympatric transplant), to np - 1 (full factorial design)
ind = 1 #ind is the number of individuals tested by transplant

sd_ind = 0	#sd_ind is the individual - level standard deviation of error
sd_pop = 10	#sd_pop is the population - level standard deviation of error
sd_hab = 8	#sd_hab is the standard deviation of habitat effects
sd_dem = 7	#sd_dem is the standard deviation of deme quality effects

SA = 10	#SA is the sympatric vs. allopatric contrast (i.e., local adaptation)
constant = 50	#constant is the baseline fitness

#Generate transplant experiment data stocked in "W" dataframe
W = make_data(np, d, ind, sd_ind, sd_pop, sd_hab, sd_dem, SA, constant)

W$habs = as.factor(W$habs)
W$dems = as.factor(W$dems)
W$symp = as.factor(W$symp)
lm_SA_pop = lm(W$fitness ~ W$habs + W$dems + W$symp)
anova(lm_SA_pop) #in this case residuals include both the rest of the G*E interaction and the error, and the default p value is correct
