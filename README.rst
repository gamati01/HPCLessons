HPCLessons
==========


1. Introduction

Repository for pdf and simple program for HPC lessons for PHD students (2025 Edition)

https:https://phd.uniroma1.it/web/OffertaFormativaErogataCiclo.aspx?c=40&i=3556&l=IT

2. Lesson structures and main arguments


* Lesson 0	(HPC-0.pdf)
	Course Presentation
		* Teachers
		* repository
		* lessons
		* pre-requisites

* Lesson 1: 01/04/2025
	introduction to HPC
	introduction to CoCalc
	Homework:
		* exploitation MM multiplication for 512^2 Matrices (blind)
	Material:
		* Books et al..
	Take home message:
		* HPC is a complicate stuff, many skills are needed
		* Good/Bad performances depends from user
		* TO exploit performance HW (vanilla) knowledge is mandatory
		
* Lesson 2: 08/04/2025
	MM multiplication verification
		* Why so different performances?
	How the HW is:
 		* Von Neumann Model (Store/Load architecture)
		* Memory Gap
	Why the memory GAP?
		* Hiding Memory gap: 
 	                How a Floating point units works
			Caching: pro/cons
	Take home message:
		* Data access is the key to go slow/fast


* Lesson 3: 15/04/25 
 	How a Floating point units works
		* Serial, Pipelined, Superpipelined, OOO
		* Vector
		* precision
		* ILP
	FP Exploition
		* blocking
		* padding
	Homework:
		* MM multiplication with FP exploitation
	Take home message:
		* Fill the FP Units its the key to performance
		* ILP
		* ---> Better Algorithm
		* ---> GPU 


* 4	29/04/2025
	MM multiplication verification
	Algorithm
		* order of complexity
	The profiling: understand where you spent more time
		* gprof
		* code instrumentation 
	Implementation: put all together
		* Data access
		* FP exploitation
		* precision
	Example
		* Crivello di Eratostene
		* TSP
		* Laplace
			* Jacobi
			* SOR
	Homework:
		* Jacobi vs. SOR. : which is the most efficente? in which region?
	Take home message
		* Algorithm chioce is a crucial step
		* some operations are (very) more expensive than others
		* mandatory to understand where you are spending more time


* 5 	05/05/2025
	Compiler 
		* what can do
		* what cannot do
		* how to force or inihbit
		* instructions vs. statement
	languages
		* C
		* Fortran
		* matlab
	Example: 	
		* reduced precision: Source of error (computing pi)
	Take home message
		* Compiler can both boost or depress performance


* 6 	12/05/2025
	CPU Parallel paradigm 1
	 	 Scalaility, Ahmdal law
		* OpenMP
			* best practices
		* MPI
			* best practices
		* intrinsic parallelism (do concurrent/....)
	Take home message
                * There are different way to exploit parallel computing
		* CPU optimization is a preliminary step before going to GPUs

* 7 	19/05/2025	
	Parallel paradigm 2 (GPU)
		* OpenACC
		* OpenMP offload
		* cuda/cuda Fortran
		* OpenCL/Sycl
	Take home message
		* GPU are important, but they are "ultima ratio regum"

* 8 	26/05/2025	
		* Conclusion/Comments
		* Performance Portability issues
		* wrap-up
		* Q/A
	Take home message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
		* HPC is th sum of different skills. Be courious

* 	Additional
		* HPC-spoiler (LBM smagorinski inplementation)

