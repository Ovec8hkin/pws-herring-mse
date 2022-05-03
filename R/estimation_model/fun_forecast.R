# fun_forecast.r
# Created by John Trochta
# Date created:  07/28/2020
# Summary:
# Run forecast from assessment to provide estimates for management procedure

fun_forecast <- function(est_nya, w_a_a, pk,
                      log_MeanAge0, annual_age0devs, sigma_age0devs, mdm_c, hydPWSSC_q, alpha_v, beta_v){
  
  # Calculate average weight over last 5 years
  
  # Calculate mean recruitment over last 10 years as the forecasted recruitment (age 3)
  
  # Project numbers-at-age using recorded catches from current year 
  # Call numbers this year using last year's info for ages 4-8
  for(int j=5;j<=nage-1;j++){
  projected_N_y_a(j)=((N_y_a(nyr_tobefit,j-1)-(SeAC(nyr_tobefit,j-1)*N_se(nyr_tobefit)+gc(nyr_tobefit,j-1)+pk*pc(nyr_tobefit,j-1)))*Sur_summer(nyr_tobefit,j-1)-fbc(nyr_tobefit,j-1))*forecast_Sur_winter(j-1);
  }

 # Calc numbers this year using last year's info for plus group
  for(int j=nage;j<=nage;j++){
    projected_N_y_a(j)=((N_y_a(nyr_tobefit,j-1)-(SeAC(nyr_tobefit,j-1)*N_se(nyr_tobefit)+gc(nyr_tobefit,j-1)+pk*pc(nyr_tobefit,j-1)))*Sur_summer(nyr_tobefit,j-1)-fbc(nyr_tobefit,j-1))*forecast_Sur_winter(j-1)+((N_y_a(nyr_tobefit,j)-(SeAC(nyr_tobefit,j)*N_se(nyr_tobefit)+gc(nyr_tobefit,j)+pk*pc(nyr_tobefit,j)))*Sur_summer(nyr_tobefit,j)-fbc(nyr_tobefit,j))*forecast_Sur_winter(j);
  }
  
  # Project early spawning biomass at age (Maturity * average weight * numbers-at-age)
  
  # Sum early SB at age to get pre-fishery biomass for input into HCR
}
  