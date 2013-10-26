CCFLAGS = -ansi -O3
LIBS = -lm
SRC_FILES = aux_fun.c buffer.c complex.c main.c notes.c onset_detection.c postponing.c spec_fun.c stdft.c
DEST = asymut

default : $(DEST)
$(DEST) : $(SRC_FILES)
	gcc $(SRC_FILES) $(CCFLAGS) $(LIBS) -o $(DEST)
clean :
	rm $(DEST) -f
