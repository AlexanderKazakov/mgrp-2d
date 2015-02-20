CC=g++
CFLAGS=-std=c++11 -c
OBJDIR=src/
SOURCES=$(OBJDIR)main.cpp $(OBJDIR)Stratum.cpp $(OBJDIR)Breaker.cpp\
$(OBJDIR)Fracture.cpp $(OBJDIR)Element.cpp $(OBJDIR)Field.cpp $(OBJDIR)util.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mhf-2d	

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ -ltinyxml -lgsl -lmgl 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:   
	rm -f $(EXECUTABLE) $(OBJECTS)