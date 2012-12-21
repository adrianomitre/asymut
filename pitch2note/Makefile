# Warning: '-march=core2' was added manually, remove if it does not apply do destination platform
#
COMP_FLAGS = -ansi -O3 -march=core2 -static
LIBS = -lm
SRC_FILES = aux_fun.c notes.c main.c
DEST = pitch2note

default : $(DEST)
$(DEST): $(SRC_FILES)
	gcc $(SRC_FILES) $(COMP_FLAGS) $(LIBS) -o $(DEST)
clean :
	rm $(DEST) -f
