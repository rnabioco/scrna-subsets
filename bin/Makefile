CC = g++
CFLAGS = -ggdb -Wall -O2 -std=c++0x 
LIBS = -lz -lm -lpthread


all: barcode_tag_bam append_read_name count_bcs rgtag_bam

.cpp:
	$(CC) $(CFLAGS) -Ihtslib -o $@ $< $(LIBS) htslib/libhts.a

