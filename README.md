# Air pollution as a driver of the epidemiological dynamics of influenza and other respiratory infections
Inter-Disciplinary Research Project --- A.Y. 2022-23

EPSRC Centre for Doctoral Training in Statistical Applied Mathematics (SAMBa)


## The Problem
While COVID has grabbed the headlines in recent years, seasonal influenza remains a significant threat and burden to global public health. There is a well-documented correlation between air pollution and influenza, along with other respiratory infections. Ozone, particulate matter and other air contaminants can inflame the lungs and compromise the immune system. Consequently, long and short-term exposure to air pollution increases the likelihood of severe disease following infection. It may also increase underlying susceptibility to viral and bacterial infections. Furthermore, some studies have found that fine particulate matter can increase virus transmission by facilitating penetration of viral particles into the lungs. See below for references. There are numerous studies of the relationship between infection incidence and air pollution in the literature. These studies generally find that poor air quality correlates with increased hospitalisation rates i.e. the incidence of severe cases. The objective of this project is to delve deeper into this relationship by developing mathematical, statistical and machine learning methods to tease apart the influences of air pollution on influenza susceptibility, transmissibility and disease severity.


## The Aim
The aim of this research is to explore the correlation between influenza cases and PM10 levels in Sweden by conducting a spatial statistical analysis. Despite being regarded as a country with relatively low levels of air pollution compared to other nations, Sweden still experiences high pollution levels in various urban areas, particularly during winter months when wood-burning stoves are commonly used for heating (IQAir 2023). Furthermore, due to transportation and industrial activities, pollutants such as nitrogen oxides and particulate matter can exceed recommended levels in some regions, posing a health hazard to certain areas and populations. On the other hand, influenza is a notable public health issue in Sweden, as epidemics of the disease are experienced every winter and result in a high incidence of sickness, hospitalizations, and even death. According to the Swedish Public Health Agency
(2023), there were over 20, 600 laboratory-confirmed influenza cases during the 2017-2018 influenza season. Hence, by investigating the spatial distribution of influenza cases and pollution levels, we can determine high-risk populations and regions, as well as examine the disease’s transmissibility in various areas.


## The Data and The Code
Real data about:

(a) pollution, from the European Environment Agency [https://www.eea.europa.eu/];

(b) weather, from the Swedish Meteorological and Hydrological Institute [https://www.smhi.se/];

(c) influenza, from the Public Health Agency of Sweden [https://www.folkhalsomyndigheten.se/].

The implementation of the influenza model is inspired by the work of *Wilson, K. & Wakefield, J. (2018), ‘Pointless spatial modeling’, Biostatistics 21(2), 17–32*, and the original code can be found in the GitHub repository [https://github.com/wilsonka/pointless-spatial-modeling].
