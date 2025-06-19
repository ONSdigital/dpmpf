# Porting R code to C++

## Introduction
This document exists as a handover of the work the author was doing, porting an R module to C++, for the purpose of estimating the population of England, Wales and Northern Ireland, using admin data: namely births, deaths, immigration and emigration. Using Bayesian stats, it should be possible to accurately estimate the population.

The C++ code that's already been ported from R exists here: https://github.com/ONSdigital/dpmaccpf/tree/main/src

Maybe you're already an experienced R programmer who also knows how to use the Rcpp tools to create high-performance C++ modules. In which case, this guide is not for you.

Maybe you're a C++ programmer - like the author - who has never done any R programming. In which case, this guide is _definitely_ for you.

If you're some way in-between the two, then this guide might prove useful, but please forgive anything which seems blindingly obvious to you, and just skim over it looking for anything interesting.

The guiding principle that the author has followed, is to develop the C++ code in an "R-like" way, so that it looks, smells, and feels familiar to R programmers, although it might not feel very "C++ like" to a pure C++ programmer. The use of Rcpp conventions is always preferred, instead of more common C++ ways of doing things, where an Rcpp option exists.

The code the author has written is as close to a like-for-like port as was possible, in the author's opinion, to better aid any debugging and QA of the code, by those reading the original R code and checking that the C++ is logically identical.

## Development Environment
All you need is R, RStudio and a C++ compiler installed.

## R and Rcpp Namespaces
There is a well-established convention that the `Rcpp` header is imported and the namespace is used, at the top of every source file. However, there also exists a `R` namespace which contains many similarly named functions, which are needed.

As a general rule, the `Rcpp` namespace contains functions which deal with vectors, and the `R` namespace contains functions which deal with scalars. The webpage the author frequently uses to check for the existence of a function is this one: https://dirk.eddelbuettel.com/code/rcpp/html/namespacemembers_func.html

## Sugar
"Sugar" is the short nickname for "syntactic sugar" meaning non-C++ conventions which work in Rcpp by means of a kind of magic: namely C++ macros, to make common "R-like" tasks very easy.

### Vectors
When we want to declare a vector in R we would write `my_numeric <- numeric(42)` and in C++ we would commonly use the `new` keyword to declare anything, but syntactic sugar allows us to write `NumericVector my_numeric (42);` which is clearly very convenient.

To access items in a vector, we can use array-like syntax sugar, but we **must remember that in C++ arrays start from element 0 not 1**. So, to store a scalar at the start of the vector, we write `my_numeric[0] = 123.456;`. Likewise we can read a value out of the vector like this: `double the_first_value = my_numeric[0];`.

We can also declare vector sequences with a start and end value, but only where the values are increasing: `Range my_positive_sequence = seq(3, 9);`. If we need an 'inverse' sequence we would have to reverse an increasing version: `IntegerVector my_inverse_sequence = rev(my_positive_sequence);`.

We can create and fill a vector with a value, like this: `NumericVector lots_of_nines = rep(9.999, 42);`

### Matrix / Matrices
To declare a matrix, we write `NumericMatrix my_matrix (2, 4);`

To access the scalar values in the matrix, we use this syntax: `double one_value = my_matrix(0, 0);`. Likewise to put a value into the matrix we do this: `my_matrix(0, 0) = 123.456;`. Note the round brackets, not square!

To fill a matrix with an initial value, the easiest way I've found that works is as follows: `std::fill(my_matrix.begin(), my_matrix.end(), 123.456);`. In my experience, the other documented way of doing it does not work.

One row of a matrix can be accessed like this: `NumericMatrix::Row my_row = my_matrix(0,_);`

Likewise one column of a matrix can be accessed like this: `NumericMatrix::Column my_column = my_matrix(_,0);`

To get the number of rows do this: `int number_of_rows = counts_true.rows();`

To get the number of columns do this: `int number_of_columns = counts_true.cols();`

## Other Rcpp Things
Rcpp has a ton of other helper constants and functions, to make it easier to do "R-like" work in C++.

### NA
To deal with NA values, there's a constant called `NA_REAL` which can be used to represent NA. To test if a value is NA or not, use the `R_IsNA()` function.

### Infinity
To deal with infinity, there are constants called `R_PosInf` and `R_NegInf`. To test if a value is infinite or not, use the `is_infinite()` function.

## Classes
Normal C++ classes can be written, but methods which are exposed to R cannot be overloaded (same method name, different parameters) because there's no way for R to infer which version of the method the caller wants to call.

The use of STL vector and other common C++ types is also strongly discouraged, in favour of the Rcpp types. For example, do not use `std::vector` but instead use `NumericVector` as a parameter, where the code is going to be called from R.

### Exposing Classes to R
Every class which is going to be used in R must be exposed by using this macro in the header file: `RCPP_EXPOSED_CLASS(MyClass)`. Further, in the source file, the class, its constructor, its methods, and any publicly accessible member variables, must also be declared with a macro block as follows:

```cpp
RCPP_MODULE(mod_myclass) {
  class_<MyClass>("MyClass")
  .constructor<int, double>("Here is a description of the class")
  .field( "myMember", &MyClass::myMember )
  .method( "myFunction", &MyClass::myFunction )
  ;
}
```

Every class is exported as a module, which must be loaded in an R file so that it can be used. An example of how to load the Rcpp module in R is: `Rcpp::loadModule("mod_myclass", TRUE)`

Once the Rcpp module is exported to R, the class can be used in R like this: `my_object <- new(MyClass, 42, 123.456)`

### Compiling Rcpp Classes
Whenever you write a new class, or modify some existing code, it's worth knowing a couple of commands to help with development.

`Rcpp::compileAttributes()` will generates the bindings required to call C++ functions from R for functions decorated with `// [[Rcpp::export]]` or otherwise exported to R from C++.

`devtools::document()` will compile the C++ without running the tests.

### Passing Lists of Class Objects to Functions
When we want to pass a list of objects (instances of classes) from R to C++, we can use the `List` type.

To get the objects out of the `List` we use the `as()` function to convert to the correct C++ type.

For example, to iterate over a list of objects and call a particular method:
```cpp
void callFunctionOnList(List myList) {
    int itemCount = myList.length();
    
    for (int i=0; i < itemCount; i++) {
      MyClass& item = as<MyClass&>(myList[i]);
      item.myFunction();
    }
}
```

### Polymorphism
Polymorphism works, but as with all object-oriented programming, you should prefer composition over inheritance. Polymorphism is generally undesirable in software engineering. You can see in the codebase that `CdmNoregBase` and `CdmWithregBase` are abstract base classes which are polymorphically implemented by subclasses. However, there are other more preferable ways this could have been implemented, so the author discourages this pattern. *JOHN BRYANT: Polymorphism is very common in R, which is why it's used here. But it would be good to discuss this, and change if necessary.*

### Memory Management / Memory Leaks
The author admits that they have no idea how the "sugar" is handling memory management, and therefore suspects that memory leaks are highly likely to occur.

Generally speaking, we should prefer that R creates all data entities, and is therefore responsible for freeing up that memory when it's no longer needed. To assist that, the convention has been, wherever possible, for C++ functions to modify data that's passed in, instead of allocating new objects and returning those.

## What's Left to do?
[John Bryant](https://github.com/johnrbryant) has a much better idea about the remaining than the author, but it seems like the Particle Filter is the large as-yet untackled area of work.

The remaining work occasionally deal with 3-dimensional matrices, which are not supported "out of the box" by Rcpp, and as such, techniques used in the [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html) package have been proposed. From a cursory glance it looks like the [Cube](http://arma.sourceforge.net/docs.html#Cube) class is a likely candidate, but the author cannot, unfortunately, offer any useful further insight.

Hopefully, a lot of the remaining porting will be a simple copy-paste exercise, using the coding conventions that have already been used extensively in the project, with the only notable exception being the aforemention 3-dimensional matrices. Mercifully, the 3D matrices appear to be only a small part of the remaining R code to be ported.

## In Conclusion
The author had no prior knowledge of R when they started porting the R code to C++ on November 27, 2021. The knowledge gleaned in the space of 46 days is captured in this document (written on January 12, 2022) and in the C++ code in the GitHub repository.

In short, this is evidence that any uninitiated programmer, with or without knowledge of R or even C++, is capable of helping this project towards completion, hopefully somewhat assisted by this document.
