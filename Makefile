CC = g++
CFLAGS = -pthread -std=c++11 -O3 -g -Wall -Wextra
INC = -Iinclude -Llib
LIBS = -lgmp -lm
RM = rm -f

DIR_OBJS=obj
DIR_SRCS=src
DIR_BINS=bin

VC_OBJS = \
$(DIR_OBJS)/field.o \
$(DIR_OBJS)/poly.o \
$(DIR_OBJS)/parameters.o \
$(DIR_OBJS)/main.o \


VC_DEPS = \
$(DIR_OBJS)/field.d \
$(DIR_OBJS)/poly.d \
$(DIR_OBJS)/parameters.d \
$(DIR_OBJS)/main.d \

VC_TARGET = $(DIR_BINS)/vc.exe


# Each subdirectory must supply rules for building sources it contributes
$(DIR_OBJS)/%.o: $(DIR_SRCS)/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CC) $(INC) $(CFLAGS) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

$(DIR_OBJS)/%.o: $(DIR_SRCS)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CC) $(INC) $(CFLAGS) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



# All Target
all: $(VC_TARGET)

# Tool invocations
$(VC_TARGET): $(VC_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC Archiver'
	$(CC) $(VC_OBJS) $(INC) $(CFLAGS) $(LIBS) -o $(VC_TARGET)
	@echo 'Finished building target: $@'
	@echo ' '


# Other Targets
clean:
	-$(RM) $(VC_OBJS) $(VC_DEPS) $(VC_TARGET)
	-@echo ' '

.PHONY: all clean dependents

