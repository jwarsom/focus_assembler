#Make file for overlapper algorithm
CC=g++
OBJECT_DIR=build
BIN_DIR=bin
CFLAGS=-c -g -Wall -Wno-char-subscripts 
LDFLAGS= 
OBJS1=$(OBJECT_DIR)/Preprocessor.o $(OBJECT_DIR)/Fragment_Index.o $(OBJECT_DIR)/Labels.o $(OBJECT_DIR)/Dictionary.o $(OBJECT_DIR)/SortFragments.o $(OBJECT_DIR)/MergeFragments.o
OBJS2=$(OBJECT_DIR)/Aligner.o $(OBJECT_DIR)/Fragment_Index.o $(OBJECT_DIR)/Seed.o $(OBJECT_DIR)/Chain.o $(OBJECT_DIR)/Dictionary.o
OBJS3=$(OBJECT_DIR)/SerialAlign.o $(OBJECT_DIR)/Fragment_Index.o $(OBJECT_DIR)/Seed.o $(OBJECT_DIR)/Chain.o $(OBJECT_DIR)/Dictionary.o

all: $(BIN_DIR)/preprocess $(BIN_DIR)/align $(BIN_DIR)/serialAlign

$(BIN_DIR)/preprocess: $(OBJS1) 
	$(CC) $(LDFLAGS) $(OBJS1) -o $(BIN_DIR)/preprocess

$(OBJECT_DIR)/Preprocessor.o: src/Preprocessor.cpp src/Fragment_Index.h src/Labels.h src/Dictionary.h src/SortFragments.h src/MergeFragments.h
	$(CC) $(CFLAGS) src/Preprocessor.cpp -o $(OBJECT_DIR)/Preprocessor.o

$(OBJECT_DIR)/Fragment_Index.o: src/Fragment_Index.cpp src/Dictionary.h
	$(CC) $(CFLAGS) src/Fragment_Index.cpp -o $(OBJECT_DIR)/Fragment_Index.o

$(OBJECT_DIR)/Labels.o: src/Labels.cpp src/Labels.h
	$(CC) $(CFLAGS) src/Labels.cpp -o $(OBJECT_DIR)/Labels.o

$(OBJECT_DIR)/SortFragments.o: src/SortFragments.cpp src/SortFragments.h
	$(CC) $(CFLAGS) src/SortFragments.cpp -o $(OBJECT_DIR)/SortFragments.o

$(OBJECT_DIR)/MergeFragments.o: src/MergeFragments.cpp src/MergeFragments.h
	$(CC) $(CFLAGS) src/MergeFragments.cpp -o $(OBJECT_DIR)/MergeFragments.o

$(BIN_DIR)/align: $(OBJS2)
	$(CC) $(LDFLAGS) $(OBJS2) -o $(BIN_DIR)/align

$(OBJECT_DIR)/Aligner.o: src/Aligner.cpp src/SuffixArray.h src/Fragment_Index.h src/MappingValues.h src/Seed.h src/Chain.h src/BandedAlignment.h src/Overlapper.h 
	$(CC) $(CFLAGS) src/Aligner.cpp -o $(OBJECT_DIR)/Aligner.o

$(OBJECT_DIR)/Seed.o: src/Seed.cpp src/Seed.h
	$(CC) $(CFLAGS) src/Seed.cpp -o $(OBJECT_DIR)/Seed.o	

$(OBJECT_DIR)/Chain.o: src/Chain.cpp src/Chain.h
	$(CC) $(CFLAGS) src/Chain.cpp -o $(OBJECT_DIR)/Chain.o	

$(BIN_DIR)/serialAlign: $(OBJS3)
	$(CC) $(LDFLAGS) $(OBJS3) -o $(BIN_DIR)/serialAlign

$(OBJECT_DIR)/SerialAlign.o: src/SerialAlign.cpp
	$(CC) $(CFLAGS) src/SerialAlign.cpp -o $(OBJECT_DIR)/SerialAlign.o

$(OBJECT_DIR)/Dictionary.o: src/Dictionary.cpp
	$(CC) $(CFLAGS) src/Dictionary.cpp -o $(OBJECT_DIR)/Dictionary.o
clean:
	rm $(OBJECT_DIR)/*.o $(BIN_DIR)/preprocess $(BIN_DIR)/align $(BIN_DIR)/serialAlign
