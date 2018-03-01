CC = g++
CFLAGS = -ggdb -Wall -O2 -std=c++0x 
LIBS = -lz -lm -lpthread


all: barcode_tag_bam append_read_name count_bcs rgtag_bam get_molecule_info

.cpp:
	$(CC) $(CFLAGS) -Ihtslib -o $@ $< bamreader.cpp $(LIBS) htslib/libhts.a

