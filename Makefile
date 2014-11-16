CC=g++
CFLAGS=-c
OBJDIR=src/
SOURCES=$(OBJDIR)main.cpp $(OBJDIR)Stratum.cpp $(OBJDIR)Fluid.cpp\
$(OBJDIR)Fracture.cpp $(OBJDIR)Break.cpp $(OBJDIR)Field.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mgrp-2d	

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ -ltinyxml -lgsl -lmgl

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:   
	rm -f $(EXECUTABLE) $(OBJECTS)