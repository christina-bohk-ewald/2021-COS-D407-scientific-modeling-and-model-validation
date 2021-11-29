#
## Basic functions for week 4
#

## 1: to interpolate IFR estimates of Verity et al. into single years of age: 

to_ungroup <- function(to_ungroup,nr_grouped_years){
	
	seq_ungrouped_years <- seq(0,length(to_ungroup)*nr_grouped_years)	
	cumsum_to_ungroup <- cumsum(c(sum(to_ungroup),to_ungroup))
	grouped_time_points <- c(0,(1:length(to_ungroup))*nr_grouped_years)
	
	applied_smooth_spline <- smooth.spline(x=grouped_time_points,y=cumsum_to_ungroup)
	predict_cumsum_ungroup <- predict(applied_smooth_spline,x=seq_ungrouped_years)$y
	ungrouped <- diff(predict_cumsum_ungroup)
	return(ungrouped)
}

to_ungroup_spar <- function(to_ungroup,nr_grouped_years,spar){
	
	seq_ungrouped_years <- seq(0,length(to_ungroup)*nr_grouped_years)	
	cumsum_to_ungroup <- cumsum(c(sum(to_ungroup),to_ungroup))
	grouped_time_points <- c(0,(1:length(to_ungroup))*nr_grouped_years)
	
	applied_smooth_spline <- smooth.spline(x=grouped_time_points,y=cumsum_to_ungroup,spar=spar)
	predict_cumsum_ungroup <- predict(applied_smooth_spline,x=seq_ungrouped_years)$y
	ungrouped <- diff(predict_cumsum_ungroup)
	return(ungrouped)
}

## 2: to ungroup remaining life expectancy:

get_ungrouped_ex_2015_2020 <- function(country_name, lt_1950_2020){
	current_period_data <- lt_1950_2020[which(lt_1950_2020[,8]=="2015-2020"),]
	current_period_data <- current_period_data[which(current_period_data[,3]==country_name),]  
	current_ex_data <- as.numeric(current_period_data[,19])
	smooth_current_ex_data <- smooth.spline(x=c(0,1,seq(5,100,5)),y=current_ex_data)
	new_x <- c(seq(0,0.99,0.01),seq(1,4.99,0.01),seq(5,100,0.01))
	predict_smooth_current_ex_data <- predict(smooth_current_ex_data,new_x,len=new_x)
	return(predict_smooth_current_ex_data)
}

## 3: to scale IFRs from RC (China, Verity et al.) onto COI: 

map_ifr_betw_ref_and_one_coi_thanatAge <- function(coi,lt_1950_2020,ungrouped_ifr_by_single_age_china_sp){

	ifr_coi_mapped_rc_china_based_on_thanat_x <- matrix(NA,nr=1,nc=length(ungrouped_ifr_by_single_age_china_sp))
	
	rownames(ifr_coi_mapped_rc_china_based_on_thanat_x) <- coi

	current_pop_insert <- coi

	for(chronAge in 1:90){
		current_ref_y <- get_ungrouped_ex_2015_2020(country_name="China",
		lt_1950_2020)$y
			
		current_ref_x <- get_ungrouped_ex_2015_2020(country_name="China", 
		lt_1950_2020)$x

		current_coi_y <- get_ungrouped_ex_2015_2020(country_name=current_pop_insert, 
		lt_1950_2020)$y
			
		current_coi_x <- get_ungrouped_ex_2015_2020(country_name=current_pop_insert, 
		lt_1950_2020)$x
				 	
		current_y_ref_of_chronAge <- current_ref_y[which(current_ref_x==(chronAge-1))]
		equal_y <- which(round(current_coi_y,3)==round(current_y_ref_of_chronAge,3))[1]
			
		if(is.na(equal_y)){
			n <- 0
			while(is.na(equal_y)){
				equal_y <- which(round(current_coi_y,3)==(round(current_y_ref_of_chronAge,3)-n))[1]
				n <- n+0.001 
			} ## while	
		} ## if
		
		equivalent_x_coi <- current_coi_x[equal_y]
		
		if((round(equivalent_x_coi,0)+1)>length(ungrouped_ifr_by_single_age_china_sp)){
			equivalent_x_coi <- 89
		}

		ifr_coi_mapped_rc_china_based_on_thanat_x[1,equivalent_x_coi] <- 
		ungrouped_ifr_by_single_age_china_sp[chronAge]

	} ## for chronAge

	return(ifr_coi_mapped_rc_china_based_on_thanat_x)

} ## function


## 3b: basic function to map ungrouped ifr based on thanatological age and ungrouped ex for France 

map_ifr_betw_assigned_ref_and_one_coi_thanatAge <- function(ref,coi,deaths,lt_1950_2020,ungrouped_ifr_by_single_age_china_sp){

	ifr_coi_mapped_rc_china_based_on_thanat_x <- matrix(NA,nr=1,nc=length(ungrouped_ifr_by_single_age_china_sp))
	
	rownames(ifr_coi_mapped_rc_china_based_on_thanat_x) <- coi

	current_pop_insert <- coi

	output_equivalent_x_coi <- c(0)

 		for(chronAge in 1:90){
			current_ref_y <- get_ungrouped_ex_2015_2020(country_name=ref, lt_1950_2020)$y
			current_ref_x <- get_ungrouped_ex_2015_2020(country_name=ref, lt_1950_2020)$x

			current_coi_y <- get_ungrouped_ex_2015_2020(country_name=current_pop_insert, lt_1950_2020)$y
			current_coi_x <- get_ungrouped_ex_2015_2020(country_name=current_pop_insert, lt_1950_2020)$x
		 	
			current_y_ref_of_chronAge <- current_ref_y[which(current_ref_x==(chronAge-1))]
			equal_y <- which(round(current_coi_y,2)==round(current_y_ref_of_chronAge,2))[1]
			
			if(is.na(equal_y)){
				n <- 0.1
				while(is.na(equal_y)){
					equal_y <- which(round(current_coi_y,2)==(round(current_y_ref_of_chronAge,2)-n))[1]
					n <- n+0.01 
				} ## while	
			} ## if

			equivalent_x_coi <- current_coi_x[equal_y]
		
			if((round(equivalent_x_coi,0)+1)>length(ungrouped_ifr_by_single_age_china_sp)){
				equivalent_x_coi <- 89
			}

			ifr_coi_mapped_rc_china_based_on_thanat_x[1,equivalent_x_coi] <- ungrouped_ifr_by_single_age_china_sp[chronAge]

			output_equivalent_x_coi[chronAge] <- equivalent_x_coi

		} ## for chronAge

	return(ifr_coi_mapped_rc_china_based_on_thanat_x)

} ## function


## 4: basic function to aggregate scaled IFRs into 10-year age groups

aggregate_mapped_ifr_10y <- function(disaggregated_mapped_ifr){
	output <- c(0)
	for(group in 1:9){
		pos <- (1+10*(group-1)):(10+10*(group-1))
		output[group] <- sum(disaggregated_mapped_ifr[pos])
	} ## for
	return(output)
} ## function

## 5: basic function to disaggregate total deaths for ONE COI into 10-year age groups based on global pattern over age 

disaggregate_deaths_one_coi_10y <- function(coi){
	deaths_by_age <- matrix(0,nr=length(seq(0,80,10)),nc=length(5:ncol(deaths)))
	rownames(deaths_by_age) <- seq(0,80,10)
	colnames(deaths_by_age) <- colnames(deaths[,5:ncol(deaths)])
	for(day in 1:length(5:ncol(deaths))){
		deaths_by_age[,day] <- deaths[which(deaths[,"Country.Region"]==coi),(day+4)] * global_age_dist_deaths$value
	} ## for
	return(deaths_by_age)
} ## function
