# BEP-SBMD-COVID19

This repository (*in the folder "notebook"*) contains the codes and logs for my **"Bachelor Eind Project" (BEP)** regarding the topic and research of COVID-19. This work will focus on the *data of COVID-19 in the Netherlands.* 


<br>


This BEP will focus on different models for the modelling and prediction of the spread
of COVID-19 in the Netherlands. Many existing epidemic models are based on the SIR
or SEIR models. Popular model for the prediction of daily cases, is to use a 3 parameters logistic
growth curve, developed by Pierre François Verhulst. However, it is also possible to use
a more dynamic model, such as the Kalman filter. Since the spread of COVID-19 in
the Netherlands can be seen as a system that is continuously changing and as a real-time
problem, using a Kalman filter might be an ideal approach.




## Software

The programming language used for this project is Python and the codes are executed on a juypter notebook. The packages used for the models can be found in the corresponding notebooks.

## Data
The data that were used for this research are included in the folder _"data"_ and are obtained from the website of **RIVM (National Institute for Public Health and the Environment)**.
The data redarding the ICU  admission is obtained from **Stichting Nice (National Intensive Care Evaluatie): https://www.stichting-nice.nl/**

## 1. SIR and SEIR model 

Both **SIR** and **SEIR** are compartemental models, which consist of differential equations, which simulates the input/output at each compartment. Which the correct estimation of the parameters in these equations, it is possible to simulate the course/spread of COVID-19. 
<br><br>
The main advantage is that these models are is independent from the actual data, thus even without proper or insufficient testing, _**the amount of people who are infected or exposed to this infectious disease, can be estimated**_. This is important, because with this information, the government, healthcare companies or hospitals can take preventive measurements (such as necessary lockdowns, making sure that there are enough medical supplies and/or that the hospitals and ICU are not overloaded).

After all, the true number of infected person can only be verified with large testing of the population. 
### SIR model 

The **SIR (Susceptible-Infected-Recovered)** model is one of the most basic compartemental model for epidemiology. This is a dynamic simulation model and can be used to simulate the course/spread (over a time range) of an infectious disease, such as COVID-19. 
<br>
<br>
As seen from the abbreviation of the model, it contains 3 different compartments: 
1. Susceptible (amount of people that can get infected)
2. Infected    ( amount of people who are alreay infected and can infected others)
3. Recovered   (amount of people recovered from the infectious disease)

Below an example of such simulation is given (in a graph):
<br>

![](/images/sir_simulation.png)

_**Note: This curve is not the actual curve that simulates the spread of COVID-19 in the Netherlands. The parameters were freely chosen for generating the plot**_

However, this model lacks complexity to simulate COVID-19. In this simulation it can be seen that eventually all "susceptibles" would be infected and later on fully recovered. Unfortunately, this is not the case with COVID-19, where there are extra compartments such as: **Exposed** and **Death**.

By adding more compartments, such as hospital admission or ICU admission, these models can be more complex and closer to "real-life problems".

### SEIR(D) model

A more releastic model, to simulate COVID-19, is achieved by adding more compartments, such as number of **Exposed** and **Deaths**. Since, not every person who are exposed to this infectious disease are immediately infected, it would be logical to add an **Exposed** compartment. Futhermore, with each infectious disease, especially the COVID-19 the death count should also be added as a compartment. 


The **SEIR (Susceptible-Exposed-Infected-Recovered)** model, does not include the death counts (although this can be seen as *recovered*, thus a modification is to add a **Death** compartment to this model, making it the **SEIR(D)** model. 

Below an example of such simulation is given (in a graph):

![](/images/seird_simulation.png)

_**Note: This curve is not the actual curve that simulates the spread of COVID-19 in the Netherlands. The parameters were freely chosen for generating the plot**_


This model can be interesting when simulating the spread of COVID-19 in a given country, because with the correctly estimated parameters, it is possible to get an estimation of the amount of people who are **Infected** or **Exposed** to COVID-19. Therefore, even with a lack of testing, it is possible to get an idea of how many people are exposed or infected. Furthermore, which such models, it is possible to see whether if the taken measurements will "flatten the curve".

## 2. Logistic Growth (3 parameters - Verhulst)

The **Logistic Growth curve** can be used to estimate the inflection point and/or predict the maximum number of cases of COVID-19. This can also be applied on the **Hospital admissions** and **ICU admissions**. 

### Approach and advantages 
This approach uses data analytics and statistics for the best fit of the curve/model. Furthermore, the actual/recorded data is used to generate this logistic growth curve. The advantage is, that it uses real data to generate and fit a model. By using real data and fitting a model, data from other countries such as China can be used to estimate to proper parameters for the logistic growth curve. This can then be used for the prediction of the number of cases of COVID-19 in the Netherlands. 

The parameters of this logistic growth curve are estimated by using **Nonlinear Least Squares Estimation**. 

The equations of this logistic growth curve is given below: 

![](/images/logistic_growth_equation.PNG)

By estimating the parameters (alpha, beta and C), whether it is using local data (the Netherlands) or from another country, it is possible to predict the **maxmimum number of cases**, thus it gives an idea when the "peak" is achieved (in how many days after the first confirmed case).

**_Note: The true maximum or the number of confirmed cases is highly dependent on the amount of tests performed_**


### Cumulative (confirmed) cases 
After estimating the parameters and fitting the model, it is possible to get an idea of the predicted maximum number of cases in the Netherlands. These models can also be used to simulate and check whether a lockdown had any effect, by fitting a model with data before the strict measurements and a model with data during these strict measurements. A difference in the **predicted maximum** can be a sign that the measurements are effective. 

The logistic growth curve and its statistics, using the data from the Netherlands, is shown below:


![](/images/log_growth_dailycases.PNG)


However, this model highly depends on the available data. Model can be tweaked by using estimated parameters obtained from data from other countries. Furthermore, it might be more interesting **the model the hospital admissions and ICU (Intensive Care Unit) admissions**. 

This is important, because it can be predicted what the **maximum intake** is and when is "peak" is achieved. With this information, it is possible to take preventive measurements (such as necessary lockdowns, making sure that there are enough medical supplies and/or that the hospitals and ICU are not overloaded).

### ~~Hospital Admissions~~ (DEPRECATED) 

**IMPORTANT**: Logistic growth is not 'suitable' for hospital & ICU admissions, because the logistic growth curve is only taking the admissions into account, but not the discharge of the patients (recovered or deceased). 
Thus, the section of **Hospital Admission** and **ICU Admissions** under Logistic Growth, **_deprecated (on Github)_**

The logistic growth curve and its statistics for hospital admissions in the Netherlands: 

![](/images/log_growth_hosp.PNG)


### ~~ICU Admissions~~ (DEPRECATED) 

The logistic growth curve and its statistics for the ICU admissions in the Netherlands: 

![](/images/log_growth_icu.PNG)

**Note: The actual data passes the predicted maximum, but these models are still valid, meaning that there exist a chance that the actual admissions can pass the predicted maximum. Therefore, it is important use these models to achieve insights when taking preventive measurements.**


### Exponential growth (3 parameters)

Although exponential growth function are quite abstract and difficult to tell when the peak is reached, by plotting the **log-scale**, a straight line is obtained. This represents the trajectory of confirmed COVID-19 cases. By plotting the actual data (in the Netherlands), it can be seen whether the actual data deviates from the trajectory. 

A downward deviation, would mean that things are getting better, but a upward deviation is a sign that more preventive measurements should be taken. 

The trajectory and its statistics, regarding the confirmed COVID-19 cases in the Netherlands: 

![](/images/exp_trajectory.PNG)