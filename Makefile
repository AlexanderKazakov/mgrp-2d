CC=g++
CFLAGS=-c
OBJDIR=src/
SOURCES=$(OBJDIR)Engine.cpp $(OBJDIR)main.cpp $(OBJDIR)Stratum.cpp $(OBJDIR)Visualization.cpp \
$(OBJDIR)Fracture.cpp $(OBJDIR)Break.cpp $(OBJDIR)Field.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mgrp-2d	

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ -ltinyxml -lgsl -lmgl -lmgl-wnd

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:   
	rm -f $(EXECUTABLE) $(OBJECTS)