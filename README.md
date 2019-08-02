# Vehicle Route Estimation Using Telematics Data
This project repository contains the R (Python code is in backup folder) implementation of the route estimation on dummy trip points dataset for learning purposes. The code showcases:

1.  Read OpenStreetMaps data, stored in a Postgres database.
2.  Read telematics trippoints data and plot them on the map, as shown in figure below.
3.  Estimate the most probable route that the drive will consider, connecting all of the trippoints

![Alt text](exmp_trippoints.jpg?raw=true? "Example Trippoints")


## Getting started
The instructions below will give you important hints on how to get the R 
functions of this project up and running on your local/remote machine, for your 
own usage, development and/or testing purposes.


### Coding golden rules
In order to collaborate to this project, please follow the following coding 
rules:
* File naming: utils_functions.R (lower case snake case),
* Function naming: CalcAverageCons (fully CamelCase),
* Variable naming: variableName (leading lower case letter followed by CamelCase),
* Constants: kConstantName (leading k, followed by CamelCase),
* Use always ``` <- ``` for value assignment and not ``` = ```,
* Use always ``` = ``` for function arguments and not ``` <- ```,
* Spacing: Place a space around all binary operators (``` = ```, ``` + ```, ``` * ```, ``` >= ```, etc.), except in function arguments,
* Never hardcode usernames, passwords and system access information in your 
developments and upload them onto the GITLab repository,
* Never upload real data onto the GITLab repository.
* Last but not least, exploit the potential of modern practices and technology...
![Alt text](out_tech.jpg?raw=true? "Outdated Technology")


### CD&CI golden rules
In order to protect code that is mature and stable (e.g. for production) from 
code that is experimental and non stable, please follow the following Continuous 
Development & Integration (CD&CI) rules:
* Always do your development on your ```LOCAL develop``` branch, not on the 
```LOCAL master``` branch,
* Check that all *Unit Tests* are passed locally before every commit,
* Add, commit and push from ```LOCAL develop``` to ```REMOTE: ORIGIN develop```, 
and never directly to ```REMOTE: ORIGIN master```,
* Once yout developments are mature and stable enough, ask a colleague for 
serving as a reviewer of your code,
* Once agreed with the reviewer, set a *merge request* from ```ORIGIN develop``` 
towards ```ORIGIN master```, nominating that reviewer by using the GITLab 
repository functionality for merging. This is done on the web browser.


### Prerequisites
In order to be able to use this project reporitory, the following R libraries 
are required:
* jsonlite,
* devtools,
* roxygen2,
* testthat,
* R.rsp,
* knitr.


### Unit testing
To perform all available *Unit Tests* of the developed functions in the project, 
run:

* the following command in an R console
```> devtools::load_all("."); print(testthat::test_dir("tests/", reporter = "summary"))```,

* or the linux shell command below
```$ Rscript -e 'devtools::load_all("."); print(testthat::test_dir("tests/", reporter = "summary"))'```,

both from the root directory of the project.

P.S: The Unit Tests are not included at the moment.


### R Package bulding
There are two types of packages that can be built out of this project, a 
*binary* package for R running on Windows OS, or a *source* package for R 
running on Linux / UNIX OS. 

* To build a *binary* package, adapt and run the following commands 
on an R console:

```> devtools::create("path/to/package/package_name")```, if not created yet

```> devtools::document(roclets=c('rd', 'collate', 'namespace'))```

```> devtools::build(binary = TRUE, args = c('--preclean'))```

* To build a *source* package, adapt and run the following commands 
on an R console:

```> devtools::create("path/to/package/package_name")```, if not created yet

```> devtools::document(roclets=c('rd', 'collate', 'namespace'))```

```> devtools::build()```


### Author, dependencies  and build version
See DESCRIPTION file.

### License
See LICENSE file.
