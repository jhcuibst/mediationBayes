# Cui, Jinhong_Work sample_2025 

## An example of R code
### This is for creating an IR function with 95% CI

```r
z = qnorm(0.975, 0, 1)
ir.func <- function(data, var1, post_cancer, cancer){
  df_numer <- data %>% group_by(!!sym(var1), !!sym(post_cancer)) %>%  summarise(n = n()) %>% ungroup()
  df_denom <- data %>% group_by(!!sym(var1)) %>% summarise(years = sum(time_to_event_year)) %>% ungroup()
  df <- df_numer %>% 
    left_join(df_denom, by = var1 ) %>%
    mutate(ir = round(100000*n/years, 1)) %>%
    mutate(r = n/years, se = sqrt((1-r)/n)) %>%
    mutate(lb = 100000*as.numeric(exp(log(r) - z*se) ),
           ub = 100000*as.numeric(exp(log(r) + z*se) )    ) %>%
    filter(!( !!sym(post_cancer) %in% c("head_neck", "Unknown", "Other"))) %>%
    mutate(rst = paste0(n, " (", ir, " [", round(lb, 1), ", ", round(ub, 1), "]", ")", sep = "")) %>%
    select(-years, -n, -r, -se, -lb, -ub) %>%
    filter(!!sym(post_cancer) == cancer) 
  return( df)
}
```

### Using the predefined function to get result
```r
ir.func(data = anal.all, var1 ="race", post_cancer =  "post_anal", cancer = "anal")
```

## Below is an example of SAS code
### Used for C=cleaning data: convert all NA to missing
```r
data stroke_free_long3; 
	set stroke_free_long2; 
	array char_vars {*} _CHARACTER_;  /* Create an array for all variables */
	   do i = 1 to dim(char_vars);    /* Loop through all variables */ 
	      if char_vars{i} = "NA" then char_vars{i} = '';  /* If the value is "NA", change it to missing */
	   end; 
	drop i; /* Drop the loop index variable */

	array num_vars {*} _NUMERIC_;  /* Create an array for all variables */
	   do j = 1 to dim(num_vars); /* Loop through all variables */ 
	      if num_vars{j} = "NA" then num_vars{j} = .;   /* If the value is "NA", change it to missing */
	   end; 
	drop j;
run;
```

### Given cleaned data, fit a mixed effect model and output results

```r
ods select ModelInfo Estimates Covtests;

proc glimmix data = stroke_free NOCLPRINT  empirical pconv=1E-5;
class id_num Gender(ref="F") REGION(ref="NONBEL");
model ace(event="1") = Age Gender Race2 visit Race2*visit REGION
			/ dist=binary link=logit solution chisq cl;* covb;
random int / type = un subject=id_num;
covtest 0;
Estimate "Interaction Race*Time" Race2*visit 1/ exp cl;
Estimate "2003-2007 Black vs White" Race2 1 / exp cl;
Estimate "2013-2016 Black vs White" Race2 1 Race2*visit 1 /exp cl;
Estimate "White 2013 vs. 2003" Race2 0 visit 1 Race2*visit 0/ exp cl;
Estimate "Black 2013 vs. 2003" Race2 0 visit 1 Race2*visit 1/ exp cl;
run;

ods select ModelInfo Estimates Covtests;
```

Author: Cui,Jinhong; jhcui@uab.edu

## An example using my R package to conduct mediaiton analysis

### Load data: 
data(example_data)

|Name | Meanings of variables in the dataset|
|-----|----------------------------|
|x    |Exposure variable.|
|mnb  |Mediator with zero-infalted negative bionomial (ZINB) distribution.|
|xm   |The interaction between exposure and mediator. User can defined it directly in the model fitting using the same way as thos in R package _glm_.|
|im   |The indicator of zero-inflated mediator, which is used for estimating decomposed Narutal Indirect Effect (NIE).|
|y   |The outcome variable.|

#### Set up priors
```r
prior.m = c(
          set_prior("student_t(3, 0, 2.5)", class = "b"),  
          set_prior("student_t(3, 0, 2.5)", class = "b", dpar = "hu") )
prior.y = set_prior("student_t(3, 0, 2.5)", class="b")
```

#### Fit Bayesian models
```r
hnb.m <- brm(bf(mnb~ x + age, hu ~ x + age), data = example_data, prior = prior.m,  
             family = hurdle_negbinomial(), chains = 4, iter = 2000)
hnb.y <- brm(y ~  x + mnb + age + im, data = example_data, prior = prior.y,  
            family = bernoulli(), chains = 4, iter = 2000)
```

#### Get outputs of mediation analysis
```r
outmed <- medbayes(hnb.m, hnb.y, mediator="mnb", treat="x", outcome = "y", ind_mediator = "im",  
                   control.value = 0, treat.value = 1)
```
