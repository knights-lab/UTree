CC=gcc
.PHONY : all

all: utree-build utree-compress utree-search, utree-build_gg utree-search_gg

utree-build: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -o BUILD PFBITS=24 IXTYPE=uint32_t

utree-compress: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -o COMPRESS PFBITS=24 IXTYPE=uint32_t

utree-search: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -o SEARCH PFBITS=24 IXTYPE=uint32_t

utree-build_gg: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -o BUILD_GG PFBITS=24 IXTYPE=uint32_t

utree-search_gg: itree.c ; $(CC) -m64 -std=gnu11 -Ofast -march=sandybridge itree.c -o SEARCH_GG PFBITS=24 IXTYPE=uint32_t

