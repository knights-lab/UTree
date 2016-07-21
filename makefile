CC=gcc
.PHONY : all

all: utree-build utree-compress utree-search, utree-build_gg utree-search_gg

utree-build: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D BUILD -D PFBITS=24 -D IXTYPE=uint32_t -o utree-build

utree-compress: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D COMPRESS -D PFBITS=24 -D IXTYPE=uint32_t -o utree-compress

utree-search: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D SEARCH -D PFBITS=24 -D IXTYPE=uint32_t -o utree-search

utree-build_gg: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D BUILD_GG -D PFBITS=24 -D IXTYPE=uint32_t -o utree-build_gg

utree-search_gg: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -D SEARCH_GG -D PFBITS=24 -D IXTYPE=uint32_t -o utree-search_gg

