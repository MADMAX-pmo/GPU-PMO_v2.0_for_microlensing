# cfitsio library path, library name
FITSPATH=-L/usr/lib64
FITSLIB=-lcfitsio
FITSINCLUDE=-I/usr/include/cfitsio
# general include path
INCLUDE=-I/usr/include 

CUDA_INSTALL_PATH = /usr/local/cuda
LDFLAGS = -L$(CUDA_INSTALL_PATH)/lib64
NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc

CULIBPATH = -lcudart -lcurand 
CC = gcc -lstdc++
#MCC = mpicc

# for GPU Tesla A100
CUFLAG = -arch=sm_80 -Xcompiler -fPIC -dc
CULIBFLAG = -arch=sm_80 -Xcompiler -fPIC -dlink

CULIB = -lcufile
CCLIB = -lcfile -lm
CFLAG = -shared -o 
CCFLAG = -fPIC -shared -o 

SRC_DIR := src
BUILD_DIR := build
PROJECT_INCLUDE := -I. -I..

CUFILES = $(addprefix $(SRC_DIR)/, calmag.cu setinterp.cu)
CFILES = $(addprefix $(SRC_DIR)/, dm.c fits_IO.c hocll.c lens_salpetermf.c)

CUOBJECTS = $(addprefix $(BUILD_DIR)/, calmag.o setinterp.o)
LINK = $(BUILD_DIR)/link.o
CULIBSO = $(BUILD_DIR)/libcufile.so
CCLIBSO = $(BUILD_DIR)/libcfile.so
EXE = $(BUILD_DIR)/fugue


all: $(BUILD_DIR) $(EXE)
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/calmag.o: $(SRC_DIR)/calmag.cu
	$(NVCC) $(CUFLAG) -c $< -o $@
$(BUILD_DIR)/setinterp.o: $(SRC_DIR)/setinterp.cu
	$(NVCC) $(CUFLAG) -c $< -o $@
$(BUILD_DIR)/link.o: $(CUOBJECTS)
	$(NVCC) $(CULIBFLAG) $^ -o $@
$(BUILD_DIR)/libcufile.so: $(CUOBJECTS) $(BUILD_DIR)/link.o
	$(CC) -shared -o $@ $(CUOBJECTS) $(BUILD_DIR)/link.o $(LDFLAGS) $(CULIBPATH)
$(BUILD_DIR)/libcfile.so: $(CFILES)
	$(CC) $^ -fPIC -shared -o $@ $(FITSINCLUDE) $(PROJECT_INCLUDE)
$(BUILD_DIR)/fugue: $(SRC_DIR)/main.c $(BUILD_DIR)/libcfile.so $(BUILD_DIR)/libcufile.so
	$(CC) $< $(PROJECT_INCLUDE) -L$(BUILD_DIR) $(CCLIB) -L$(BUILD_DIR) $(CULIB) -o $@ $(FITSPATH) $(FITSLIB)

#all:
#	$(NVCC) $(CUFLAG) $(CUFILES)
#	$(NVCC) $(CULIBFLAG) $(CUOBJECTS) -o $(LINK)
#	$(CC) $(CFLAG) $(CULIBSO) $(CUOBJECTS) $(LINK) $(LDFLAGS) $(CULIBPATH)
#	$(CC) $(CFILES) $(CCFLAG) $(CCLIBSO) $(FITSINCLUDE) $(PROJECT_INCLUDE)
#	$(CC) $(SRC_DIR)/main.c $(PROJECT_INCLUDE) -L. $(CCLIB) -L. $(CULIB)  -o $(EXE) $(FITSPATH) $(FITSLIB)
clean:
	rm -rf $(BUILD_DIR)
	rm -f nohup.out
	
