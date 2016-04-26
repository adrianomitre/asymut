CCFLAGS = -ansi -O3
LIBS = -lm
SRC_FILES = aux_fun.c notes.c main.c
DEST = pitch2note

default : $(DEST)
$(DEST): $(SRC_FILES)
	gcc $(SRC_FILES) $(CCFLAGS) $(LIBS) -o $(DEST)
clean :
	rm $(DEST) -f
