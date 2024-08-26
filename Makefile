## 0. compile object
TARGET := $(notdir $(shell pwd))

## 1. compile configure
FFLAGS = # 编译选项 
LDFLAGS = # 链接的库位置

ifeq ($(COMPIFLE),)
   COMPIFLE = intel
endif

ifeq ($(COMPIFLE), gnu)
   AR = ar
   FC = gfortran
   # 增强警告
   FFLAGS += -Wall -Werror -Wall -Warray-bounds
   # 良好的编译习惯
   FFLAGS += -Wshadow -Wconversion
   FFLAGS += -Wunreachable-code -Wunused-parameter -Wunused-variable
endif

ifeq ($(COMPIFLE), intel)
   AR = ar
   FC = ifort
   # Enable Intel compiler warnings
   FFLAGS += -warn all
   # Good coding practices
   FFLAGS +=
   # Additional warnings
   FFLAGS += 
endif

ifeq ($(MAKECMDGOALS), $(filter $(MAKECMDGOALS),gdb vscode))
    FFLAGS += -g
else
    FFLAGS += -Ofast # Optimization flags
endif

ifneq ($(MAKECMDGOALS), $(filter $(MAKECMDGOALS), clean docs))
  $(info # Building $(TARGET) [$(FC) $(FFLAGS)])
endif

## 2. directory
# Directory Structure(FPM)
WORK_DIR := $(shell pwd)
SRC_DIR  := $(WORK_DIR)/src
APP_DIR  := $(WORK_DIR)/app
DST_DIR  := $(WORK_DIR)/build/make
INC_DIR  := $(WORK_DIR)/build/make

# Create the destination directory (`./build`)
$(shell mkdir -p $(DST_DIR))

ifeq ($(COMPIFLE), gnu)
   FFLAGS += -J$(DST_DIR)
endif

ifeq ($(COMPIFLE), intel)
   FFLAGS += -module $(DST_DIR)
endif

FFLAGS += -I$(INC_DIR)

# execuble file
EXE = $(WORK_DIR)/build/$(TARGET)
LIB = $(WORK_DIR)/build/lib$(TARGET).a

# Source files
SRCS := $(wildcard $(SRC_DIR)/*.[fF]90) # 返回匹配到的文件列表
APPS := $(wildcard $(APP_DIR)/*.[fF]90) # 返回匹配到的文件列表

# Object files
NAMES := $(notdir $(SRCS))
LIB_OBJS := $(addprefix $(DST_DIR)/, $(NAMES:%=%.o)) # 变量值的替换
NAMES := $(notdir $(APPS))
APP_OBJS := $(addprefix $(DST_DIR)/, $(NAMES:%=%.o)) # 变量值的替换

## 3. compile
$(shell ./utils/deps $(SRCS) > $(DST_DIR)/Makefile.lib.dep)
include $(DST_DIR)/Makefile.lib.dep

# link the program
$(EXE): $(LIB) $(APP_OBJS)
	@echo "link $(EXE)"
	@$(FC) $(FFLAGS) $(LDFLAGS) $(APP_OBJS) $(LIB) -o $@

# 静态模式替换
$(APP_OBJS): $(DST_DIR)/%.o: $(APP_DIR)/%
	@echo FC $@
	@$(FC) $(FFLAGS) -c $< -o $@

$(LIB): $(LIB_OBJS)
	@echo "link $(LIB)"
	@$(AR) -sr $(LIB) $^

# 静态模式替换
$(LIB_OBJS): $(DST_DIR)/%.o: $(SRC_DIR)/%
	@echo FC $@
	@$(FC) -fPIC $(FFLAGS) -c $< -o $@

## 4. run this program 
run: $(EXE)
	$(EXE)

test: $(EXE)
	$(EXE)

# Debug the program with GDB
gdb: $(EXE)
	gdb $(EXE)

# Debug the program with GDB in vscode
vscode: $(EXE)

# docs
docs:
	ford $(TARGET).md -o docs

# Clean up
clean:
	rm -rf $(DST_DIR) $(EXE) $(LIB)

.PHONY: docs clean
