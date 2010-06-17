RM := rm -rf

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
./main.cpp \
./node.cpp \
./tree.cpp \
./tree_reader.cpp \
./tree_utils.cpp 

OBJS += \
./main.o \
./node.o \
./tree.o \
./tree_reader.o \
./tree_utils.o 

CPP_DEPS += \
./main.d \
./node.d \
./tree.d \
./tree_reader.d \
./tree_utils.d 

OPT_FLAGS += \
-O3 -ffast-math -ftree-vectorize -g3 -Wall
#-g -O0

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ $(OPT_FLAGS) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# All Target
all: phyx_rt phyx_pl

# Tool invocations
phyx_rt: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -o "phyx_rt" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

phyx_pl: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -o "phyx_pl" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) phyx
	-@echo ' '

