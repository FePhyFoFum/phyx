RM := rm -rf

#this is really nothing but for testing the capabilities of the copper library
TEST_SRC = main.cpp
TEST_OBJ = main.o
TEST_DEP = main.d

#this should be the basic copper library sources, objects and deps
CPP_SRCS += \
./node.cpp \
./tree.cpp \
./tree_reader.cpp \
./tree_utils.cpp 

CPP_OBJS += \
./node.o \
./tree.o \
./tree_reader.o \
./tree_utils.o 

CPP_DEPS += \
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

# test target
test: phyx_test

phyx_test: $(TEST_OBJ) $(CPP_OBJS) 
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -o "phyx_test" $(TEST_OBJ) $(CPP_OBJS) $(CPP_LIBS)
	@echo ' '

# Tool invocations
phyx_rt: $(CPP_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -o "phyx_rt" $(CPP_OBJS) $(CPP_LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

phyx_pl: $(CPP_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -o "phyx_pl" $(CPP_OBJS) $(CPP_LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CPP_OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) phyx
	-@echo ' '

