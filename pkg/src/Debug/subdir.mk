################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../RnewMat.cpp \
../bayesMfp.cpp \
../combinatorics.cpp \
../dataStructure.cpp \
../hyperg.cpp \
../sort.cpp 

C_SRCS += \
../const.c \
../hyp2f1.c \
../mtherr.c 

OBJS += \
./RnewMat.o \
./bayesMfp.o \
./combinatorics.o \
./const.o \
./dataStructure.o \
./hyp2f1.o \
./hyperg.o \
./mtherr.o \
./sort.o 

C_DEPS += \
./const.d \
./hyp2f1.d \
./mtherr.d 

CPP_DEPS += \
./RnewMat.d \
./bayesMfp.d \
./combinatorics.d \
./dataStructure.d \
./hyperg.d \
./sort.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


