# Warning: '-march=core2' was added manually, remove if it does not apply do destination platform
#
COMP_FLAGS = -ansi -O3 -march=core2 -static
LIBS = -lm
SRC_FILES = aux_fun.c buffer.c complex.c main.c notes.c onset_detection.c postponing.c spec_fun.c stdft.c
DEST = asymut

default : $(DEST)
$(DEST) : $(SRC_FILES)
	gcc $(SRC_FILES) $(COMP_FLAGS) $(LIBS) -o $(DEST)
clean :
	rm $(DEST) -f
