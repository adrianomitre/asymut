CCFLAGS = -ansi -O3
LIBS = -lm
SRC_FILES = aux_fun.c filters.c main.c
DEST = peak_filter

default : $(DEST)
$(DEST): $(SRC_FILES)
	gcc $(SRC_FILES) $(CCFLAGS) $(LIBS) -o $(DEST)
clean :
	rm $(DEST) -f
