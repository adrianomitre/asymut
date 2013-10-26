CCFLAGS = -ansi -O3
LIBS = -lm
SRC_FILES = aux_fun.c f0_estimation_methods.c main.c
DEST = pda

default : $(DEST)
$(DEST) : $(SRC_FILES)
	gcc $(SRC_FILES) $(CCFLAGS) $(LIBS) -o $(DEST)
clean :
	rm $(DEST) -f
