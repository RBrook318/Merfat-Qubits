# Compiler and flags
F90 = gfortran
FFLAGS = -g
FLINKS = -llapack -lblas
VPATH = ../src

# Directories
SRC_DIR = src
BUILD_DIR = build
MOD_DIR = mod

# Target executable
TARGET = su2.exe

# Source files
SRC_FILES = Merfat_arr_prm.f90 Merfat_zqp_fullbasis.f90 \
            Merfat_zqp_0c_state.f90 Merfat_zqp_1c_state.f90 \
            Merfat_input.f90 Merfat_overlap.f90 Merfat_overlap_zinit_zi.f90 Merfat_init_conditions.f90 \
            Merfat_CCF.f90 Merfat_hamiltonian.f90 Merfat_drivs.f90 Merfat_step_t.f90 \
            Merfat_norm.f90 Merfat_ACF.f90 Merfat_calculation.f90  \
            Merfat_main_file.f90

# Object files
OBJ_FILES = $(addprefix $(BUILD_DIR)/,$(SRC_FILES:.f90=.o))

# Default target
all: $(TARGET)

# Link step
$(TARGET): $(OBJ_FILES)
	$(F90) $(FFLAGS) -o $@ $^ $(FLINKS)

# Compile steps
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR) $(MOD_DIR)
	$(F90) $(FFLAGS) -c $< -o $@ -J$(MOD_DIR)

# Ensure build and mod directories exist
$(BUILD_DIR) $(MOD_DIR):
	mkdir -p $@

# Clean target
clean:
	rm -rf $(BUILD_DIR) $(MOD_DIR) *.exe *~ ../run/*.exe ../src/*~ ../run/*~

# .PHONY targets
.PHONY: all clean
