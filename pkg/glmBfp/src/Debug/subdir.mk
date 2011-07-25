################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../bfgs.o \
../combinatorics.o \
../dataStructure.o \
../design.o \
../distributions.o \
../evalZdensity.o \
../fpUcHandling.o \
../functionWraps.o \
../glmBayesMfp.o \
../iwls.o \
../laplace.o \
../linalgInterface.o \
../optimize.o \
../sampleGlm.o \
../zdensity.o 

CPP_SRCS += \
../bfgs.cpp \
../combinatorics.cpp \
../dataStructure.cpp \
../design.cpp \
../distributions.cpp \
../evalZdensity.cpp \
../fpUcHandling.cpp \
../functionWraps.cpp \
../glmBayesMfp.cpp \
../iwls.cpp \
../linalgInterface.cpp \
../optimize.cpp \
../sampleGlm.cpp \
../zdensity.cpp 

OBJS += \
./bfgs.o \
./combinatorics.o \
./dataStructure.o \
./design.o \
./distributions.o \
./evalZdensity.o \
./fpUcHandling.o \
./functionWraps.o \
./glmBayesMfp.o \
./iwls.o \
./linalgInterface.o \
./optimize.o \
./sampleGlm.o \
./zdensity.o 

CPP_DEPS += \
./bfgs.d \
./combinatorics.d \
./dataStructure.d \
./design.d \
./distributions.d \
./evalZdensity.d \
./fpUcHandling.d \
./functionWraps.d \
./glmBayesMfp.d \
./iwls.d \
./linalgInterface.d \
./optimize.d \
./sampleGlm.d \
./zdensity.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/lib/R/include -I"/home/daniel/R/forge/glmBfp/src" -I/home/daniel/R/library/Rcpp/include -I/home/daniel/R/library/RcppArmadillo/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


