CC = g++
CFLAGS = -ggdb -Wall -O2 -std=c++0x 
LIBS = -lbz2 -lz -lm -lpthread


all: barcode_tag_bam append_read_name count_bcs

.cpp:
	$(CC) $(CFLAGS) -Ihtslib -o $@ $< $(LIBS) htslib/libhts.a

