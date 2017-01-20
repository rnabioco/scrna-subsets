CC = g++
CFLAGS = -ggdb -Wall -O2 -std=c++0x 
LIBS = -lbz2 -lz -lm -lpthread

barcode_tag_bam: barcode_tag_bam.cpp
	$(CC) $(CFLAGS) -Ihtslib -o $@ $< $(LIBS) htslib/libhts.a

append_read_name: append_read_name.cpp
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

dplex_illumina: dplex_illumina.cpp
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
